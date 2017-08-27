#' Differential Methylation with Classical Target Parameters
#'
#' Computes the Targeted Minimum Loss-Based Estimate of the Average Treatment
#' Effect (ATE) or the Risk Ratio, treating DNA methylation as an outcome (Y)
#' and the indicated variable of interest (which ought to be binarized) as a
#' treatment/exposure (A), using the neighbors of a given CpG site as the
#' adjustment set (W).
#'
#' @param target_site Numeric ...
#' @param methytmle_screened An object of class \code{methytmle}...
#' @param var_of_interest ...
#' @param type Character ...
#' @param corr Numeric ...
#' @param obs_per_covar Numeric ...
#' @param target_param Character ...'
#' @param g_lib Character or vector of characters...
#' @param Q_lib Character or vector of characters...
#' @param family Character ...
#' @param return_ic Logical ...
#'
#' @importFrom minfi getBeta getM
#' @importFrom tmle tmle
#' @importFrom stats cor
#'
#' @export
#'
methyvim_tmle <- function(target_site,
                          methytmle_screened,
                          var_of_interest,
                          type = c("Beta", "Mval"),
                          corr,
                          obs_per_covar,
                          target_param = c("ate", "rr"),
                          g_lib = c("SL.mean", "SL.glm"),
                          Q_lib = c("SL.mean", "SL.glm"),
                          family = c("gaussian", "binomial"),
                          return_ic = FALSE
                         ) {
  ### check arguments where possible
  #type <- match.arg(type)
  #family <- match.arg(family)

  ### get neighboring site
  in_cluster <- which(methytmle_screened@clusters %in%
                      methytmle_screened@clusters[target_site])

  ### remove target site from the set of neighbors
  only_neighbors <- in_cluster[in_cluster != target_site]

  # get expression measures based on input
  if (type == "Beta") {
    expr <- minfi::getBeta(methytmle_screened)
  } else if (type == "Mval") {
    expr <- minfi::getM(methytmle_screened)
  }

  # get measures at the target site
  y <- as.numeric(as.matrix(expr[target_site, , drop = FALSE]))

  # perform scaling of outcome if using binomial error family
  if (family == "binomial") {
    a <- min(y, na.rm = TRUE)
    b <- max(y, na.rm = TRUE)
    y_star <- (y - a) / (b - a)
  } else {
    y_star <- y
  }

  ### are there enough neighbors for this estimate to be meaningful
  if (length(only_neighbors) != 0) {
    # how many neighbors were there originally?
    n_neighbors_total <- length(only_neighbors)

    # extract measures for neighboring sites
    w <- as.data.frame(expr[only_neighbors, , drop = FALSE])

    # quick sanity check of the size of w
    stopifnot(nrow(w) == length(only_neighbors))

    # find maximum number of covariates that can be in W
    w_max <- round(length(y_star) / obs_per_covar)

    # remove neighbors that are highly correlated with target site
    if (sum(abs(stats::cor(y, t(w))) > corr) > 0) {
      # find neighbors that are highly correlated with the target site
      neighbors_corr_high <- which(abs(stats::cor(y, t(w))) > corr)

      # if all neighbors are too highly correlated, we'll simply ignore W
      if (length(neighbors_corr_high) == length(only_neighbors)) {
        w_no_corr <- NULL
        w_in <- as.data.frame(t(rep(1, length(y))))
      } else if (length(neighbors_corr_high) != 0) {
        w_no_corr <- TRUE
        w_in <- as.data.frame(w[-neighbors_corr_high, ])
      }
    } else {
      w_no_corr <- TRUE
      w_in <- w
    }

    # use PAM to reduce W by selecting medoids
    if (!is.null(w_no_corr) & nrow(w_in) > w_max) {
      message("PAM will be used to reduce W but is not yet implemented.")
      # write utility function to perform PAM clustering and select medoids
      # TODO: w <- cluster_w_pam(w)
    }

    # strictly enforces the assumption of positivity by discretizing W
    if (!is.null(w_no_corr)) {
      w_pos <- force_positivity(var_of_interest, t(w_in), pos_min = 0.15)
    } else {
      w_pos <- as.data.frame(t(w_in))
    }

    # get length of remaining neighbors in the adjustment set
    if (!is.null(w_no_corr)) {
      n_neighbors_reduced <- ncol(w_pos)
    } else {
      n_neighbors_reduced <- 0
    }

    # maximum correlation among neighbors in the adjustment set
    max_corr_w <- max(cor(y, t(w)))

    # compute the ATE
    out <- tmle::tmle(Y = as.numeric(y_star),
                      A = as.numeric(var_of_interest),
                      W = as.data.frame(w_pos),
                      Q.SL.library = Q_lib,
                      g.SL.library = g_lib,
                      family = family,
                      verbose = FALSE
                     )
  } else {
    n_neighbors_total <- 0
    n_neighbors_reduced <- 0
    max_corr_w <- NA

    # perform estimation with W = 1 if there are no neighbors
    w_int <- as.data.frame(rep(1, length(var_of_interest)))

    # compute the ATE
    out <- tmle::tmle(Y = as.numeric(y_star),
                      A = as.numeric(var_of_interest),
                      W = as.data.frame(w_int),
                      Q.SL.library = Q_lib,
                      g.SL.library = g_lib,
                      family = family,
                      verbose = FALSE
                     )
  }

  # get the influence curve estimates if so requested
  if (return_ic) {
    ate_ic <- out$estimates$IC$IC.ATE
    ate_g <- out$g$g1W
    ate_Q <- out$Qstar
    ic <- list(ate_ic, ate_g, ate_Q)
  }

  # extract and rescale estimates
  if (target_param == "ate") {
    est <- out$estimates$ATE
    est_raw <- c(est$CI[1], est$psi, est$CI[2], est$var.psi, est$pvalue)
    est_rescaled <- est_raw[1:3] * (b - a)
    var_rescaled <- est_raw[4] * ((b - a)^2)
    res <- c(est_rescaled, var_rescaled, est_raw[5], n_neighbors_total,
             n_neighbors_reduced, max_corr_w)
  } else if (target_param == "rr" & family == "binomial") {
    est <- out$estimates$RR
    est_ci_log <- log(est$CI)
    est_out <- c(est_ci_log[1], est$log.psi, est_ci_log[2], est$var.log.psi,
                 est$pvalue)
    res <- c(est_out, n_neighbors_total, n_neighbors_reduced, max_corr_w)
  } else if (target_param == "rr" & family == "gaussian") {
    stop("The Relative Risk is not estimable with a gaussian error family.")
  }
  return(res)
}

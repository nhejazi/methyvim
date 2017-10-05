#' Differential Methylation with Classical Target Parameters
#'
#' Computes the Targeted Minimum Loss-Based Estimate of the Average Treatment
#' Effect (ATE) or the Risk Ratio, treating DNA methylation as an outcome (Y)
#' and the indicated variable of interest (which ought to be binarized) as a
#' treatment/exposure (A), using the neighbors of a given CpG site as the
#' adjustment set (W). INTERNAL USE ONLY.
#'
#' @param target_site Numeric indicating the column containing the screened
#'        CpG site indices that will be looped over in TMLE procedure.
#' @param methytmle_screened An object of class \code{methytmle} with clustered
#'        sites based on genomic windows containing the screened CpG site
#'        indices column.
#' @param var_of_interest Numeric indicating the column index of the binarized
#'        variable of interest, treated as an exposure.
#' @param type Character indicating the particular measure of DNA methylation to
#'        be used as the observed data in the estimation procedure, either Beta
#'        values or M-values. The data are accessed via \code{minfi::getBeta} or
#'        \code{minfi::getM}.
#' @param corr Numeric indicating the maximum correlation that a neighboring
#'        site can have with the target site.
#' @param obs_per_covar Numeric indicating the number of observations needed for
#'        for covariate included in W for downstream analysis. This ensures the
#'        data is sufficient to control for the covariates.
#' @param target_param Character indicating the target causal parameter for
#'        which an estimator will be constructed and computed via targeted
#'        minimum loss-based estimation. Currently, this is limited to the
#'        Average Treatment Effect (ATE) and the Risk Ratio (RR), with routines
#'        from the \code{tmle} package being used for the computation.
#' @param g_lib Character or vector of characters indicating the algorithms to
#'        be implemented in SuperLearner if \code{tmle_type} is set to "glm".
#' @param Q_lib Character or vector of characters indicating the algorithms to
#'        be implemented in SuperLearner if \code{tmle_type} is set to "sl".
#' @param family Character indicating the distribution to be implemented to
#'        describe the error distribution for regressions, generally "gaussian"
#'        for a continuous outcome and "binomial" for a binary outcome.
#' @param return_ic Logical indicating whether an influence curve estimate
#'        should be returned for each site that passed through the filter.
#'
#' @return An \code{data.frame} containing the results of the Targeted Minimum
#'         Loss-based Estimation (TMLE) procedure for the target parameter of
#'         interest for a single CpG site, computed via the \code{tmle} package.
#'
#' @keywords internal
#'
#' @importFrom minfi getBeta getM
#' @importFrom stats cor
#' @importFrom cluster pam
#' @importFrom tmle tmle
#
methyvim_tmle <- function(target_site,
                          methytmle_screened,
                          var_of_interest,
                          type = c("Beta", "Mval"),
                          corr,
                          obs_per_covar,
                          target_param = c("ate", "rr"),
                          g_lib = c("SL.mean", "SL.glm", "SL.glm.interaction"),
                          Q_lib = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                          family = c("gaussian", "binomial"),
                          return_ic = FALSE
                         ) {
  ### check arguments where possible
  type <- match.arg(type)
  target_param <- match.arg(target_param)
  family <- match.arg(family)

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
      w_pam <- cluster::pam(x = w_in, k = w_max, diss = FALSE)
      w_in <- as.data.frame(w_pam$medoids)
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
    #ate_ic <- out$estimates$IC$IC.ATE
    #ate_g <- out$g$g1W
    #ate_Q <- out$Qstar
    #ic <- list(ate_ic, ate_g, ate_Q)
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
  }
  return(res)
}


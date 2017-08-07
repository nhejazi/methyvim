#' Differential Methylation with the Average Treatment Effect
#'
#' Computes the Targeted Minimum Loss-Based Estimate of the Average Treatment
#' Effect (ATE), treating DNA methylation as an outcome (Y) and the indicated
#' variable of interest (which ought to be binarized) as a treatment/exposure
#' (A), using the neighbors of a given CpG site as the adjustment set (W).
#'
#' @param methy_tmle An object of class \code{methytmle}...
#'        ...
#' @param catch_inputs_ate List...
#'        ...
#'

methyvim_ate <- function(methy_tmle,
                         catch_inputs_ate) {

  # make sure that the outcome data is of class numeric
  var_of_interest <- as.numeric(colData(methy_tmle)[, catch_inputs$var])
  if (class(var_of_interest) != "numeric") {
    var_of_interest <- as.numeric(var_of_interest)
  }

  # get names of sites to be added to output object
  cpg_screened_names <- names(methy_tmle[methy_tmle@screen_ind])

  # object of screened CpG site indices to loop over in TMLE procedure
  methy_tmle_ind <- seq_along(methy_tmle_screened@screen_ind)

  #=============================================================================
  # perform TMLE estimation of L-ATE for all genomic (CpG) sites individually
  # ============================================================================
  methy_vim_out <- foreach::foreach(i_site = methy_tmle_ind,
                                    .packages = c("tmle", "class", "gtools"),
                                    .combine = rbind) %dopar% {

    ### create message for monitoring
    message(paste("Beginning estimation for site", i_site))

    ### get target site
    target_site <- methy_tmle_screened@screen_ind[i_site]

    ### get neighboring site
    in_cluster <- which(methy_tmle_screened@clusters %in%
                        methy_tmle_screened@clusters[target_site])

    ### remove target site from the set of neighbors
    only_neighbors <- in_cluster[in_cluster != target_site]

    # get expression measures based on input
    if (catch_inputs_ate$type == "Beta") {
      expr <- minfi::getBeta(methy_tmle_screened)
    } else if (catch_inputs_ate$type == "Mval") {
      expr <- minfi::getM(methy_tmle_screened)
    }

    # get measures at the target site and perform scaling
    y <- as.numeric(as.matrix(expr[target_site, , drop = FALSE]))
    a <- min(y, na.rm = TRUE)
    b <- max(y, na.rm = TRUE)
    y_star <- (y - a) / (b - a)

    ### are there enough neighbors for this estimate to be meaningful
    if (length(only_neighbors) != 0) {
      # how many neighbors were there originally?
      n_neighbors_total <- length(only_neighbors)

      # extract measures for neighboring sites
      w <- as.data.frame(expr[only_neighbors, , drop = FALSE])

      # quick sanity check of the size of w
      stopifnot(nrow(w) == length(only_neighbors))

      # find maximum number of covariates that can be in W
      w_max <- round(length(y_star) / catch_inputs_ate$obs_per_covar)

      # remove neighbors that are highly correlated with target site
      if (sum(abs(cor(y, t(w))) > catch_inputs_ate$corr) > 0) {
        # find neighbors that are highly correlated with the target site
        neighbors_corr_high <- which(abs(cor(y, t(w))) > catch_inputs_ate$corr)

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
      out <- tmle(Y = as.numeric(y_star),
                  A = as.numeric(var_of_interest),
                  W = as.data.frame(w_pos),
                  Q.SL.library = catch_inputs_ate$tmle_args$Q_lib,
                  g.SL.library = catch_inputs_ate$tmle_args$g_lib,
                  family = catch_inputs_ate$tmle_args$family,
                  verbose = FALSE
                 )
    } else {
      n_neighbors_total <- 0
      n_neighbors_reduced <- 0
      max_corr_w <- NA

      # perform estimation with W = 1 if there are no neighbors
      w_int <- as.data.frame(rep(1, length(var_of_interest)))

      # compute the ATE
      out <- tmle(Y = as.numeric(y_star),
                  A = as.numeric(var_of_interest),
                  W = as.data.frame(w_int),
                  Q.SL.library = catch_inputs_ate$tmle_args$Q_lib,
                  g.SL.library = catch_inputs_ate$tmle_args$g_lib,
                  family = catch_inputs_ate$tmle_args$family,
                  verbose = FALSE
                 )
    }

    # get the influence curve estimates if so requested
    if (catch_inputs_ate$return_ic) {
      ate_ic <- out$estimates$IC$IC.ATE
      ate_g <- out$g$g1W
      ate_Q <- out$Qstar
      ic <- list(ate_ic, ate_g, ate_Q)
    }

    # extract and rescale estimates
    est <- out$estimates$ATE
    est_raw <- c(est$CI[1], est$psi, est$CI[2], est$var.psi, est$pvalue)
    est_rescaled <- est_raw[1:3] * (b - a)
    var_rescaled <- est_raw[4] * ((b - a)^2)
    res <- c(est_rescaled, var_rescaled, est_raw[5], n_neighbors_total,
             n_neighbors_reduced, max_corr_w)
  }

  #=============================================================================
  # TMLE procedure is now done, so let's just make the output object pretty...
  # ============================================================================
  methy_vim_out <- as.data.frame(methy_vim_out)
  colnames(methy_vim_out) <- c("lower_CI_ATE", "est_ATE", "upper_CI_ATE", "Var",
                               "pval", "n_neighbors_all", "n_neighbors_w",
                               "max_corr_w")
  rownames(methy_vim_out) <- cpg_screened_names
  return(methy_vim_out)
}

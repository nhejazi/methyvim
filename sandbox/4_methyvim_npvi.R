# screen CpG sites using LIMMA method
methy_tmle_screened <- limma_screen(methytmle = methy_tmle,
                                    var_int = catch_inputs_npvi$var,
                                    type = catch_inputs_npvi$type,
                                    cutoff = catch_inputs_npvi$filter_cutoff)

# NOTE: work out NPVI procedure "by hand"
# create clusters
methy_tmle_screened <- cluster_sites(methy_tmle = methy_tmle_screened,
                                     window_size = catch_inputs_npvi$window)

# extract the variable of interest and ensure numeric
var_of_interest <- colData(methy_tmle_screened)[, catch_inputs_npvi$var]
if (class(var_of_interest) != "numeric") {
  var_of_interest <- as.numeric(var_of_interest)
}

# setup NPVI params
confid_level <- catch_inputs_npvi$tmle_args$npvi_conf
descr <- catch_inputs_npvi$tmle_args$npvi_descr
alpha <- (1 - confid_level)
samp_size <- length(var_of_interest)

# run the NPVI procedure
methy_tmle_ind <- seq_along(methy_tmle_screened@screen_ind)
sites_rand <- sort(sample(methy_tmle_ind, 50))
sites <- names(methy_tmle_screened[methy_tmle_screened@screen_ind[sites_rand],])

methy_vim_out <- foreach::foreach(i_site = sites_rand,
                                  .packages = c("tmle.npvi", "SuperLearner"),
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
  if (catch_inputs_npvi$type == "Beta") {
    expr <- minfi::getBeta(methy_tmle_screened)
  } else if (catch_inputs_npvi$type == "Mval") {
    expr <- minfi::getM(methy_tmle_screened)
  }

  ### are there enough neighbors for this estimate to be meaningful
  if (length(only_neighbors) != 0) {
    # how many neighbors were there originally?
    n_neighbors_total <- length(only_neighbors)

    # get measures at the target site and cutoff for null value in NPVI
    x <- as.numeric(as.matrix(expr[target_site, ]))
    cutoff <- quantile(abs(x), probs = catch_inputs_npvi$tmle_args$npvi_cutoff)
    tx_zero <- which(abs(x) < cutoff)

    # extract measures for neighboring sites
    w <- as.data.frame(expr[only_neighbors, , drop = FALSE])

    # quick sanity check of the size of w
    stopifnot(nrow(w) == length(only_neighbors))

    # find maximum number of covariates that can be in W
    w_max <- round(length(var_of_interest) / catch_inputs_npvi$obs_per_covar)

    # remove neighbors that are highly correlated with target site
    if (sum(abs(cor(x, t(w))) > catch_inputs_npvi$corr) > 0) {
      # find neighbors that are highly correlated with the target site
      neighbors_corr_high <- which(abs(cor(x, t(w))) > catch_inputs_npvi$corr)

      # if all neighbors are too highly correlated, we'll simply ignore W
      if (length(neighbors_corr_high) == length(only_neighbors)) {
        w_no_corr <- NULL
        w_in <- as.data.frame(t(rep(1, length(x))))
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

    # maximum correlation among neighbors in the adjustment set
    max_corr_w <- max(cor(x, t(w)))

    # set the values below the quantile cutoff to zero
    x[tx_zero] <- 0

    # strictly enforces the assumption of positivity by discretizing W
    #if (!is.null(w_no_corr)) {
    #  x_bin <- as.numeric(x != 0)
    #  w_pos <- force_positivity(x_bin, t(w_in), pos_min = 0.15)
    #} else {
    #  w_pos <- as.data.frame(t(w_in))
    #}

    # get length of remaining neighbors in the adjustment set
    if (!is.null(w_no_corr)) {
      n_neighbors_reduced <- nrow(w_in)
    } else {
      n_neighbors_reduced <- 0
    }

    # create observed data matrix for input into tmle.npvi
    obs_data_in <- as.data.frame(cbind(var_of_interest, x, t(w_in)))
    colnames(obs_data_in) <- c("Y", "X", paste0("W", 1:(ncol(obs_data_in) - 2)))

    # compute the NPVI
    out <- tmle.npvi(obs = obs_data_in,
                     f = descr$f,
                     flavor = "superLearning",
                     stoppingCriteria = descr$stoppingCriteria,
                     cvControl = descr$cvControl,
                     nMax = descr$nMax,
                    )
  } else {
    n_neighbors_total <- 0
    n_neighbors_reduced <- 0
    max_corr_w <- NA

    # get measures at the target site and cutoff for null value in NPVI
    x <- as.numeric(as.matrix(expr[target_site, ]))
    cutoff <- quantile(abs(x), probs = catch_inputs_npvi$tmle_args$npvi_cutoff)
    tx_zero <- which(abs(x) < cutoff)
    x[tx_zero] <- 0

    # perform estimation with W = 1 if there are no neighbors
    # NOTE: apparently need at least 2 W columns for NPVI to function
    w_int <- as.data.frame(replicate(2, rep(1, length(var_of_interest))))

    # create observed data matrix for input into tmle.npvi
    obs_data_in <- as.data.frame(cbind(var_of_interest, x, w_int))
    colnames(obs_data_in) <- c("Y", "X", paste0("W", 1:(ncol(obs_data_in) - 2)))

    # compute the NPVI
    out <- tmle.npvi(obs = obs_data_in,
                     f = descr$f,
                     flavor = "superLearning",
                     stoppingCriteria = descr$stoppingCriteria,
                     cvControl = descr$cvControl,
                     nMax = descr$nMax,
                    )
  }

  # get the influence curve estimates if so requested
  if (catch_inputs_npvi$return_ic) {
    npvi_ic <- out$estimates$IC$IC.ATE
    npvi_g <- out$g$g1W
    npvi_Q <- out$Qstar
    ic <- list(npvi_ic, npvi_g, npvi_Q)
  }

  # change confidence level if not default of 95%
  if (confid_level != 0.95) {
    setConfLevel(out, confid_level)
  }

  # extract results from NPVI for a SINGLE SITE
  # NOTE: it's unclear how to properly extract the variance from NPVI
  npvi_est <- getPsi(out)
  pval <- getPValue(out)
  names(pval) <- NULL
  CI = npvi_est + c(-1,1) * getSic(out) * qnorm(1 - alpha/2) / sqrt(samp_size)
  res <- c(CI[1], npvi_est, CI[2], pval)
  out <- c(res, n_neighbors_total, n_neighbors_w, max_corr_w)
}
methy_vim_out <- as.data.frame(methy_vim_out)
colnames(methy_vim_out) <- c("lower_CI_ATE", "est_ATE", "upper_CI_ATE", "pval",
                             "n_neighbors_all", "n_neighbors_w", "max_corr_w")
rownames(methy_vim_out) <- sites

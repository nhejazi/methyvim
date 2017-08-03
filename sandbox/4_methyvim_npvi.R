# screen CpG sites using LIMMA method
methy_tmle_screened <- limma_screen(methytmle = methy_tmle,
                                    var_int = catch_inputs$var,
                                    type = catch_inputs$type,
                                    cutoff = catch_inputs$cutoff)

# NOTE: work out NPVI procedure "by hand"
## first, randomly generate an outcome since we don't have this in the data
outcome_int <- rnorm(nrow(colData(methy_tmle)))

# create clusters
methy_tmle_screened <- cluster_sites(methy_tmle = methy_tmle_screened)

# find all cases that have no missing values
cases_complete <- complete.cases(colData(methy_tmle_screened))

# remove subjects with missingness from the outcome
if (length(outcome_int) > sum(cases_complete)) {
  y <- as.numeric(outcome_int[cases_complete])
} else {
  y <- as.numeric(outcome_int)
}

# run the NPVI procedure
confid_level = 0.95
descr <- catch_inputs$npvi_descr

# setup NPVI params
alpha <- (1 - confid_level)
samp_size <- length(y)

# setup shortened TMLE loop
methy_tmle_ind <- seq_along(methy_tmle_screened@screen_ind)
sites_rand <- sort(sample(methy_tmle_ind, 50))

methy_vim_out <- foreach::foreach(i_site = sites_rand,
                                  .packages = c("tmle.npvi", "SuperLearner"),
                                  .combine = rbind) %dopar% {

  # get target site
  target_site <- methy_tmle_screened@screen_ind[i_site]

  # get neighboring site
  in_cluster <- which(methy_tmle_screened@clusters %in%
                      methy_tmle_screened@clusters[target_site])

  # remove target site from the set of neighbors
  only_neighbors <- in_cluster[in_cluster != target_site]

  # are there enough neighbors for this estimate to be meaningful
  if (length(only_neighbors) != 0) {
    # get expression measures based on input
    if (catch_inputs$type == "Beta") {
      expr <- minfi::getBeta(methy_tmle_screened)
    } else if (catch_inputs$type == "Mval") {
      expr <- minfi::getM(methy_tmle_screened)
    }

    # get measures at the target site and cutoff for null value in NPVI
    x <- expr[target_site, ]
    cutoff <- quantile(x, probs = catch_inputs$npvi_cut)

    # extract measures for neighboring sites
    w <- expr[only_neighbors, ]

    # find maximum number of covariates that can be in W
    w_max <- round(length(x) / catch_inputs$obs_per_var)

    # use HOPACH to force dimension reduction of W
    if (nrow(w) > w_max) {
      if (catch_inputs$dm) {
        message("HOPACH has not yet been implemented for dimension reduction.")
      } else {
        # select the least correlated neighbors for adjustment set
        w <- w[order(abs(cor(y, t(w))))[seq_len(w_max)], ]

        # remove neighbors that are highly correlated if still too many
        if (nrow(w) > w_max) {
          if (sum(abs(cor(x, t(w))) > catch_inputs$corr) > 0) {
            w <- w[-which(abs(cor(x, t(w))) > catch_inputs$corr), ]
          }
        }
      }
    }

    # set values below cutoff to zero for NPVI
    x[which(x < cutoff)] <- 0

    # create observed data matrix for input into tmle.npvi
    obs_data_in <- as.data.frame(cbind(y, x, t(w)))
    colnames(obs_data_in) <- c("Y", "X", paste0("W", 1:(ncol(obs_data_in) - 2)))

    # compute the ATE
    out <- tmle.npvi(obs = obs_data_in,
                     f = descr$f,
                     flavor = "superLearning",
                     stoppingCriteria = descr$stoppingCriteria,
                     cvControl = descr$cvControl,
                     nMax = descr$nMax,
                    )

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
  } else {
    message(paste("No valid neighbors found in window around site", i_site))
    no_neighbors <- c(no_neighbors, i_site)
    res <- rep(NA, 4)
  }
}

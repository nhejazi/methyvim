# screen CpG sites using LIMMA method
methy_tmle_screened <- limma_screen(methytmle = methy_tmle,
                                    var_int = catch_inputs$var,
                                    type = catch_inputs$type)

# NOTE: work out ATE procedure "by hand"
var_of_interest <- as.numeric(colData(methy_tmle_screened)[, catch_inputs$var])
if (class(var_of_interest) != "numeric") {
  var_of_interest <- as.numeric(var_of_interest)
}

# create clusters
methy_tmle_screened <- cluster_sites(methy_tmle = methy_tmle_screened)

# find all cases that have no missing values
cases_complete <- complete.cases(colData(methy_tmle_screened))

# run the ATE procedure
no_neighbors <- NULL
methy_tmle_ind <- seq_along(methy_tmle_screened@screen_ind)
methy_vim_out <- foreach::foreach(i_site = methy_tmle_ind,
                                  .packages = c("tmle", "class", "gtools"),
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

    # get measures at the target site and perform scaling
    y <- expr[target_site, ]
    a <- min(y, na.rm = TRUE)
    b <- max(y, na.rm = TRUE)
    y_star <- (y - a) / (b - a)

    # extract measures for neighboring sites
    w <- expr[only_neighbors, ]

    # find maximum number of covariates that can be in W
    w_max <- round(length(y_star) / catch_inputs$subj_per_covar)

    # remove neighbors that are highly correlated with target site
    w <- w[-which(abs(cor(y, t(w))) > catch_inputs$corr_max), ]

    # use HOPACH to force dimension reduction of W
    if (nrow(w) > w_max) {
      if (catch_inputs$dm) {
        message("HOPACH has not yet been implemented for dimension reduction.")
      } else {
        # select the least correlated neighbors for adjustment set
        w <- w[order(abs(cor(y, t(w))))[seq_len(w_max)], ]
      }
    }

    # strictly enforces the assumption of positivity by discretizing W
    w_pos <- force_positivity(var_of_interest, t(w))

    # compute the ATE
    out <- tmle(Y = as.numeric(y_star),
                A = as.numeric(var_of_interest),
                W = as.data.frame(w_pos),
                Q.SL.library = catch_inputs$Q_lib,
                g.SL.library = catch_inputs$g_lib,
                family = catch_inputs$family,
                verbose = FALSE
               )
    # extract and rescale estimates
    est <- out$estimates$ATE
    est_raw <- c(est$CI[1], est$psi, est$CI[2], est$var.psi, est$pvalue)
    est_rescaled <- est_raw[1:3] * (b - a)
    var_rescaled <- est_raw[4] * ((b - a)^2)
    res <- c(est_rescaled, var_rescaled, est_raw[5])
  } else {
    message(paste("No valid neighbors found in window around site", i_site))
    no_neighbors <- c(no_neighbors, i_site)
    res <- rep(NA, 5)
  }
}
colnames(methy_vim_out) <- c("lower_ci", "ATE", "upper_ci", "Var", "pval")

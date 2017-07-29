#' Differential Methylation Analysis with Nonparametric Variable Importance
#'
#' Computes the Targeted Minimum Loss-Based Estimate of...
#'
#' @param data_grs - ...
#'        ...
#' @param var_int - ...
#'        ...
#' @param min_sites Numeric indicating the minimum number of sites to be kept
#'        after the application of a filtering procedure.
#' @param neighbors Numeric indicating the size (in bp) of genomic windows to
#'        use when assigning CpG sites to groups.
#' @param parallel c(logical, numeric) - whether to use or number of cores to be
#'        used in parallelizing the TMLE-based procedure across genomic sites.
#' @param family (character) - specification of error family; available options
#'        are "binomial" or "gaussian" (default).
#' @param g_lib - library of learning algorithms to be used in fitting the model
#'        to "g" in the Targeted Minimum Loss-Based Estimation procedure.
#' @param Q_lib - library of learning algorithms to be used in fitting the model
#'        to "Q" in the Targeted Minimum Loss-Based Estimation procedure.
#'
#' @export methyvim
#'

methyvim <- function(data_grs,
                     var_int = 1,
                     cpg_is = "exposure",
                     type = c("Beta", "M"),
                     vim = c("ate", "npvi"),
                     neighbors = 1e3,
                     normalize = NULL,
                     filter = TRUE,
                     min_sites = 1e4,
                     family = "gaussian",
                     g_lib = c("SL.mean", "SL.glm", "SL.randomForest"),
                     Q_lib = c("SL.mean", "SL.randomForest"),
                     parallel = TRUE,
                     return_ic = TRUE,
                     shrink_ic = TRUE
                     ) {

  # ============================================================================
  # catch input and return in output object for user convenience
  # ============================================================================
  call <- match.call(expand.dots = TRUE)
  type <- match.arg(type)

  # ============================================================================
  # catch inputs to pass to downstream functions (for estimation and such)
  # ============================================================================
  catch_inputs <- list(data = data_grs, var = var_int, cpg_iss = cpg_is,
                       type = type, vim = vim, neighbors = neighbors,
                       normalize = normalize, filter = filter,
                       min_sites = min_sites, family = family,
                       g_lib = g_lib, Q_lib = Q_lib, parallel = parallel,
                       return_ic = return_ic, shrink_ic = shrink_ic)
  catch_inputs <- check_inputs(catch_inputs)

  # ============================================================================
  # invoke S4 class constructor for "methadapt" object
  # ============================================================================
  methy_tmle <- .methytmle(catch_inputs[["data"]])
  methy_tmle@call <- call

  #=============================================================================
  # set up parallelization if so desired
  # ============================================================================
  if (catch_inputs[["parallel"]]) {
    set_parallel()
  }

  # ============================================================================
  # operate on the type of data specified
  # ============================================================================
  if (catch_inputs[["type"]] == "Beta") {
    cpg_data <- getBeta(methy_tmle)
  } else if (catch_inputs[["type"]] == "M") {
    cpg_data <- getM(methy_tmle)
  }

  # ============================================================================
  # normalize array data if so requested
  # ============================================================================
  if (!is.null(catch_inputs[["normalize"]])) {
    #cpg_data <- ...
  }
  if (catch_inputs[["type"]] == "Beta") {
    getBeta(methy_tmle) <- cpg_data
  } else if (catch_inputs[["type"]] == "M") {
    getM(methy_tmle) <- cpg_data
  }

  #=============================================================================
  # check if there is missing data in the phenotype-level matrix and drop if so
  # ============================================================================
  cols_na <- colSums(sapply(colData(methy_tmle), is.na))
  na_var_id <- names(which(cols_na != 0))
  message(paste("Missing data detected: Dropping variable(s)", na_var_id,
                "from phenotype matrix."))
  colData(methy_tmle)[na_var_id] <- NULL

  #=============================================================================
  # TMLE procedure for targeted differential methylation analysis
  # ============================================================================
  methyTMLEout <- foreach(site = 1:length(sitesReduced),
                          .combine = cbind) %dopar% {

    cluster <- as.numeric(sitesClust[site, 1, drop = FALSE])
    targetSite <- as.numeric(sitesClust[site, -1, drop = FALSE])
    nearbySites <- subset(sitesClust, sitesClust[, 1] == cluster)
    nearbySites <- as.data.frame(t(nearbySites[, -1, drop = FALSE]))

    if(ncol(nearbySites) > round(nrow(nearbySites) / perSubj)) {
      print(paste("insufficient sample size to compute L-ATE for site", site))
      out <- rbind(Inf, Inf, Inf, Inf, Inf)
    } else {
      print(paste("Estimating L-ATE for", site, "of", length(sitesReduced)))
      out <- tmle(Y = targetSite,
                  A = exposureDiscrete,
                  W = nearbySites,
                  Q.SL.library = Q_lib,
                  g.SL.library = g_lib,
                  family = "gaussian",
                  verbose = FALSE
                 )

      out <- out$estimates$ATE
      out <- rbind(out$CI[1], out$psi, out$CI[2], out$var.psi, out$pvalue)
    }
  }
  methytmle <- as.data.frame(t(methyTMLEout))
  colnames(methytmle) <- c("lowerCI", "ATE", "upperCI", "Var", "pval")
  methadapt$tmleOut <- methytmle
  return(methadapt)
}

#' Differential Methylation Analysis via Targeted Minimum Loss-Based Estimation
#'
#' Computes the Targeted Minimum Loss-Based Estimate of the Local Average
#' Treatment Effect (L-ATE) for each of a reduced set of genomic sites.
#'
#' @param GRanges - ...
#'        ...
#' @param varEffect - ...
#'        ...
#' @param sitesReduced - ...
#'        ...
#' @param parallel (logical, numeric) - whether to use or number of cores to be
#'        used in parallelizing the TMLE-based procedure across genomic sites.
#' @param family (character) - specification of error family; available options
#'        are "binomial" or "gaussian" (default).
#' @param Q_lib - library of learning algorithms to be used in fitting the model
#'        to "Q" in the Targeted Minimum Loss-Based Estimation procedure.
#' @param g_lib - library of learning algorithms to be used in fitting the model
#'        to "g" in the Targeted Minimum Loss-Based Estimation procedure.
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @export biomarkertmle
#'

methytmle <- function(GRanges,
                      varEffect,
                      sitesReduced,
                      parallel = TRUE,
                      family = "gaussian",
                      Q_lib = c("SL.mean", "SL.randomForest"),
                      g_lib = c("SL.mean", "SL.glm", "SL.randomForest")
                     ) {

  # ============================================================================
  # catch input and return in output object for user convenience
  # ============================================================================
  call <- match.call(expand.dots = TRUE)

  # ============================================================================
  # invoke S3 class constructor for "biotmle" object
  # ============================================================================
  methadapt <- methadapt(call = call, tmleOut = NULL)

  #=============================================================================
  # set up parallelization based on input
  # ============================================================================
  if (class(parallel) == "numeric") doParallel::registerDoParallel(parallel)
  if (class(parallel) == "logical") {
    nCores <- parallel::detectCores()
    if (nCores > 1) {
      doParallel::registerDoParallel(nCores)
    } else {
      warning("option 'parallel' is set to TRUE but only 1 core detected.")
    }
    if (parallel == FALSE) {
      warning("parallelization has been set to FALSE: the estimation procedure
               will likely take on the order of days to run to completion.")
    }
  }

  #=============================================================================
  # organize GRanges input object based on reduced sites to be considered
  # ============================================================================
  # keep only the complete cases and normalize...
  sitesCompleteCases <- as.data.frame(mcols(cpg_gr))[, keep_cases[[1]]]
  sitesClust <- as.data.frame(cbind(clusters, sitesCompleteCases))

  #=============================================================================
  # discretize effect variable of interest to define target causal parameter
  # ============================================================================
  exposure <- designMats_min[[1]][, 2]
  exposureDiscrete <- as.numeric(exposure > quantile(exposure)["25%"])

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

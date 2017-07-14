methytmle_disc <- function(sumExp,
                           clusters,
                           outcomeVar,
                           targetSites,
                           txCutoff = "50%",
                           type = "exposure",
                           numberSites = NULL,
                           positivityMin = 10,
                           covarPerSubj = 20,
                           parallel = TRUE,
                           family = "binomial",
                           g_lib = list("SL.mean", "SL.glm", "SL.bayesglm"),
                           Q_lib = list("SL.mean", "SL.glm", "SL.knn", "SL.gam")
                          ) {

  # ============================================================================
  # check inputs are of the correct type before proceeding
  # ============================================================================
  type <- match.arg(type)
  family <- match.arg(family)

  # ============================================================================
  # catch input and return in output object for user convenience
  # ============================================================================
  call <- match.call(expand.dots = TRUE)

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
        warning("parallelization set to FALSE: manually abort procedure.")
     }
  }

  #=============================================================================
  # heuristics for getting the value of the outcome variable
  # ============================================================================
  if (length(outcomeVar) == 1) {
    # assume referring to column of design matrix if a scalar
    outcomeValues <- as.data.frame(colData(sumExp))[, outcomeVar]
  } else if (length(outcomeVar) > 1) {
    # assume actual outcome data if a vector
    outcomeValues <- outcomeVar
  } else {
    warning("inappropriate value specified for 'outcomeVar' argument")
  }

  # make sure that the outcome data is of class numeric
  if (class(outcomeValues) != "numeric") {
    warning("outcome not numeric...coercing to numeric, but may cause errors.")
    outcomeValues <- as.numeric(outcomeValues)
  }

  # find all cases that have no missing values
  cases_complete <- complete.cases(colData(sumExp))

  # scale outcome variable for use with binomial TMLE
  a = min(outcomeValues, na.rm = TRUE)
  b = max(outcomeValues, na.rm = TRUE)
  y_star <- (outcomeValues - a) / (b - a)

  # remove all missing values in case outcome specified via design matrix
  if(length(y_star) > sum(cases_complete)) {
    y_star <- as.numeric(y_star[cases_complete])
  }

  #=============================================================================
  # set number of sites to be tested...
  # ============================================================================
  if (is.null(numberSites)) {
    numberSites <- length(targetSites)
  }

  #=============================================================================
  # perform TMLE estimation of L-ATE for all genomic (CpG) sites individually
  # ============================================================================
  methyTMLEout <- foreach::foreach(site = 1:numberSites,
                                   .packages = c("tmle", "class", "gtools"),
                                   .combine = rbind) %dopar% {
     set.seed(64014617)
     target <- targetSites[site]
     cluster <- as.numeric(clusters[target])

     # create vector for target site + matrix for all others WITH cluster ID
     targetSite <- as.numeric(as.data.frame(assay(sumExp)[target, ]))
     methData <- as.data.frame(cbind(clusters, as.data.frame(assay(sumExp))))

     # find neighbors based on clusters and remove target itself from neighbors
     nearbySites <- subset(methData, methData[, 1] == cluster)
     nearbySites <- as.data.frame(t(nearbySites[, -1, drop = FALSE]))
     targetIndex <- as.numeric(which(colSums(nearbySites - targetSite) == 0))
     nearbySites <- as.data.frame(nearbySites[, -targetIndex, drop = FALSE])

     ## rarely, the only neighbor a target site has is itself. In these cases,
     ## it is obviously not possible to define the TMLE for the parameter of
     ## interest (since W = NULL).
     if (dim(nearbySites)[2] == 0) {
       res <- rep(NA, 7)
     } else {
       # remove subjects for which data is not "complete" (based on design)
       targetSite <- targetSite[cases_complete]
       nearbySites <- as.data.frame(nearbySites[cases_complete, ])

       # discretize gene of interest for use in Local ATE
       targetScaled <- as.numeric(scale(targetSite))
       cutoffTargetScaled <- quantile(abs(targetScaled))[as.character(txCutoff)]
       targetInd <- as.numeric(abs(targetScaled) > cutoffTargetScaled)

       # discretize measures from nearby sites, dropping bad sites if necessary
       posMinProp = round(positivityMin / length(targetInd), digits = 3)
       niceNearbySites <- optimPosit(A = targetInd,
                                     W = nearbySites,
                                     posMin = posMinProp)

       if (dim(niceNearbySites)[2] == 0) {
         resNA <- rep(NA, 5)
         numNeighbors <- 0
         res <- c(numNeighbors, cutoffTargetScaled, resNA)
       } else {
         # find maximum number of neighboring sites that can be controlled for
         # (based on a heuristic centered around sample size)
         maxNearbySites <- round(nrow(niceNearbySites) / covarPerSubj)

         # TMLE procedure to estimate variable importance of CpG sites
         if (ncol(niceNearbySites) < maxNearbySites) {
           print(paste("Computing Local ATE for site", site, "/", numberSites))
           out <- tmle(Y = as.numeric(y_star),
                       A = as.numeric(targetInd),
                       W = as.data.frame(niceNearbySites),
                       Q.SL.library = Q_lib,
                       g.SL.library = g_lib,
                       family = family,
                       verbose = FALSE
                      )
           numNeighbors <- ncol(as.data.frame(niceNearbySites))
         } else {
           print(paste("Computing Local ATE for site", site, "/", numberSites))
           maxNiceNearbySites <- niceNearbySites[, 1:maxNearbySites]
           out <- tmle(Y = as.numeric(y_star),
                       A = as.numeric(targetInd),
                       W = as.data.frame(maxNiceNearbySites),
                       Q.SL.library = Q_lib,
                       g.SL.library = g_lib,
                       family = family,
                       verbose = FALSE
                      )
           numNeighbors <- ncol(as.data.frame(maxNiceNearbySites))
         }
         est <- out$estimates$ATE
         est_raw <- c(est$CI[1], est$psi, est$CI[2], est$var.psi, est$pvalue)
         est_rescaled <- est_raw[1:3] * (b - a)
         var_rescaled <- est_raw[4] * ((b - a)^2)
         res <- c(numNeighbors, cutoffTargetScaled, est_rescaled, var_rescaled,
                  est_raw[5])
       }
     }
  }
  methTMLE <- as.data.frame(methyTMLEout)
  colnames(methTMLE) <- c("Wj", "TxCutoff", "lowerCI", "LocalATE", "upperCI",
                          "Variance", "pvalue")
  rownames(methTMLE) <- names(rowRanges(sumExp)[targetSites[1:numberSites]])

  # use the FDR-MSA adjustment procedure from Tuglus & van der Laan (2008)
  totalTests <- nrow(assay(sumExp))
  tmle_pvals <- methTMLE$pvalue
  tmle_pvals[which(is.na(tmle_pvals))] <- 1
  resultsFDR <- FDR_msa(pvals = tmle_pvals, totalTests = totalTests)
  methTMLE$pvalFDR <- resultsFDR

  # build SummarizedExperiment object for output
  tmleOut_rR <- subset(rowRanges(sumExp), names(sumExp) %in% rownames(methTMLE))
  tmleOut_cD <- colData(sumExp)[complete.cases(colData(sumExp)), ]
  tmleOut_se <- SummarizedExperiment(rowRanges=tmleOut_rR, colData=tmleOut_cD)
  rowData(tmleOut_se) <- methTMLE
  metadata(tmleOut_se) <- list(type = "methadapt output", created = Sys.time(),
                               call = call)
  # return output as SummarizedExperiment container object
  return(tmleOut_se)
}

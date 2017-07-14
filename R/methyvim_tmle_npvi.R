methytmle_cont <- function(sumExp,
                           clusters,
                           outcomeVar,
                           targetSites,
                           type = "exposure",
                           confid_level = 0.95,
                           nullMode = "support",
                           nullSupport = 0.20,
                           nullRange = c(-2, 2),
                           numberSites = NULL,
                           parallel = TRUE,
                           cvControl = 2,
                           nMax = 30,
                           ...
                          ) {

  # ============================================================================
  # check inputs are of the correct type before proceeding
  # ============================================================================
  type <- match.arg(type)

  # ============================================================================
  # catch input and return in output object for user convenience
  # ============================================================================
  call <- match.call(expand.dots = TRUE)

  # ============================================================================
  # create TMLE-NPVI description list for convenient input later
  # ============================================================================
  descr <- list(f = identity, iter = 10, cvControl = cvControl, nMax = nMax,
                stoppingCriteria = list(mic = 0.001, div = 0.001, psi = 0.01))

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
    # NOTE: add r2weight here?
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

  # remove all missing values if necessary
  if (length(outcomeValues) > sum(cases_complete)) {
    y <- as.numeric(outcomeValues[cases_complete])
  } else {
    y <- as.numeric(outcomeValues)
  }


  # ============================================================================
  # get confidence level and sample size
  # ============================================================================
  alpha <- (1 - confid_level)
  sampSize <- length(y)

  #=============================================================================
  # set number of sites to be tested...(just pick the first n sites)
  # ============================================================================
  if (is.null(numberSites)) {
    numberSites <- length(targetSites)
  }

  #=============================================================================
  # perform TMLE estimation of L-ATE for all genomic (CpG) sites individually
  # ============================================================================
  methyTMLEout <- foreach::foreach(site = 1:numberSites,
                                   .packages = c("tmle.npvi", "SuperLearner"),
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
     ## it is obviously not possible to define a TMLE for the parameter of
     ## interest (since W = NULL).
     if (dim(nearbySites)[2] == 0) {
       res <- rep(NA, 6)
     } else {
       # remove subjects for which data is not "complete" (based on design)
       targetSite <- targetSite[cases_complete]
       nearbySites <- as.data.frame(nearbySites[cases_complete, ])

       # set values in the null range to 0 for the target site
       if (nullMode == "scientific") {
         nullCount <- sum(targetSite > nullRange[1] & targetSite < nullRange[2])
         if (nullCount < round(0.2 * length(targetSite))) {
           warning("fewer than 20% of target sites in the null range.")
           message(paste("NPVI estimation procedure unstable for site", site))
           badSite <- TRUE
         } else {
           badSite <- FALSE
         }
         nullObs <- which(targetSite > nullRange[1] & targetSite < nullRange[2])
         targetSite[nullObs] <- 0
       } else if (nullMode == "support") {
         nullObs <- which(abs(targetSite) < quantile(abs(targetSite),
                                                         nullSupport))
         targetSite[nullObs] <- 0
         badSite <- FALSE
       }
     }

       # set up CpG-specific matrix for use with TMLE-NPVI
       #y <- (y - max(y)) / (max(y) - min(y))
       siteDataIn <- as.data.frame(cbind(y, targetSite, nearbySites))
       colnames(siteDataIn) <- c("Y", "X",
                                 paste0("W", 1:(ncol(siteDataIn) - 2)))

       # run TMLE-NPVI using the nearly created input data structure
       npviOut <- tmle.npvi(obs = siteDataIn,
                            f = descr$f,
                            flavor = "superLearning",
                            stoppingCriteria = descr$stoppingCriteria,
                            cvControl = descr$cvControl,
                            nMax = descr$nMax,
                           )

       # change confidence level if not default of 95%
       if (confid_level != 0.95) {
         setConfLevel(npviOut, confid_level)
       }

       # extract results from "npviOut" for a SINGLE SITE
       paramEst <- getPsi(npviOut)
       paramVar <- (getPsiSd(npviOut)^2)
       pval <- getPValue(npviOut)
       names(pval) <- NULL
       CI = paramEst + c(-1,1)*getSic(npviOut)*qnorm(1-alpha/2)/sqrt(sampSize)
       res <- c(badSite, CI[1], paramEst, CI[2], paramVar, pval)
     }
  }
  methTMLE <- as.data.frame(methyTMLEout)
  colnames(methTMLE) <- c("stableNPVI", "lowerCI", "estNPVI", "upperCI",
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

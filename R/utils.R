#' Clustering of sites to generate CpG neighborhoods
#'
#' @importFrom bumphunter clusterMaker
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom IRanges ranges
#'

cluster_sites <- function(methy_tmle, window_size = 1000) {
  gr <- SummarizedExperiment::rowRanges(methy_tmle)
  pos <- BiocGenerics::start(IRanges::ranges(gr))
  clusters <- bumphunter::clusterMaker(chr = GenomeInfoDb::seqnames(gr),
                                       pos = pos,
                                       assumeSorted = FALSE,
                                       maxGap = window_size)
  methy_tmle@clusters <- as.numeric(clusters)
  return(methy_tmle)
}

################################################################################

#' FDR-MSA correction
#'
#' Modified FDR Controlling Procedure for Multi-Stage Analyses (MJ van der Laan
#' and C Tuglus, 2009, <doi:10.2202/1544-6115.1397>)
#'
#' @importFrom stats p.adjust
#'

fdr_msa <- function(pvals, total_obs) {
  pvals_not_tested <- rep(1, total_obs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fdr_adj <- p.adjust(pvals_all, method = "fdr")
  fdr_out <- fdr_adj[seq_len(pvals)]
  return(fdr_out)
}

################################################################################

#' Set up parallelization
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'

set_parallel <- function(parallel) {
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
}

################################################################################

#' Find correlations between neighboring CpGs
#'
#' Determine those CpG sites highly correlated with a particular site, towards
#' the purpose of removing such sites from the baseline covariates included in
#' TMLEs of a parameter of interest, to avoid practical positivity violations.
#'
#'
#'
#'
#'

corr_cpg <- function() {
  message("This method has not been implemented yet.")
}

################################################################################

#' Discretize a vector
#'
#' Discretizes a non-factor input vector and returns the result as numeric.
#'
#' @param x A vector containing arbitrary data.
#'
#' @return A numeric vector with the data re-coded to based on the quantiles.
#'
#' @importFrom gtools quantcut
#'
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' discrete_by_quantile(x)

discrete_by_quantile <- function(x) {
  if( class(x) != "factor" ) {
    as.numeric(gtools::quantcut(x))
  } else {
    as.numeric(x)
  }
}

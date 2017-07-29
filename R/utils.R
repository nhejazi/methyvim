#' Clustering of sites to generate CpG neighborhoods
#'
#' @importFrom bumphunter clusterMaker
#'

cluster_sites <- function(granges, ...) {
  if(dim(mcols(granges))[2] != 0) {
    mcols(granges) <- NULL
  }
  clusters <- bumphunter::clusterMaker(chr = seqnames(granges),
                                       pos = start(ranges(granges)),
                                       assumeSorted = FALSE,
                                       maxGap = 1000)
  return(clusters)
}

################################################################################

#' FDR-MSA correction (van der Laan and Tuglus, 2009, <doi:>)
#'
#' @importFrom stats p.adjust
#'

fdr_msa <- function(pvals, total_obs) {
  pvals_not_tested <- rep(1, total_obs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fdr_adj <- p.adjust(pvals_all, method = "fdr")
  fdr_out <- fdr_adj[1:length(pvals)]
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

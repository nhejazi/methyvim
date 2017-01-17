# utility functions not for direct use by the user

#' clustering of sites for TMLE-based evaluation of neighboring regions
#'
#' @importFrom bumphunter boundedClusterMaker
#'
cluster_sites <- function(granges, ...) {
  if(dim(mcols(granges))[2] != 0) {
    mcols(granges) <- NULL
  }
  clusters <- bumphunter::boundedClusterMaker(chr = seqnames(granges),
                                              pos = start(ranges(granges)),
                                              assumeSorted = FALSE,
                                              maxClusterWidth = 3000,
                                              maxGap = 300)
  return(clusters)
}

#' FDR-MSA correction (van der Laan and Tuglus)
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

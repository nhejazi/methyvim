#' CpG Neighborhoods from Genomic Distance
#'
#' Clustering of CpG sites to define CpG neighborhoods based on distance (bp).
#' INTERNAL USE ONLY.
#'
#' @param methytmle Object of class \code{methytmle} produced from an object of
#'  class \code{GenomicRatioSet}, but containing extra slots.
#' @param window_size Numeric giving the number of base pairs used to define
#'  neighborhoods along the genome. CpG sites within this distance (bp) of one
#'  another are denoted as neighbors. Chromosome boundaries and other biological
#'  constraints are respected (see the documentation available for
#'  \code{bumphunter::clusterMaker} for details).
#'
#' @return An object of class \code{methytmle} with the "clusters" slot filled
#'  in. The "clusters" slot contains a \code{numeric} vector as long as the
#'  number of CpG sites. Each entry in the vector is a neighborhood assignment
#'  used in the estimation procedure.
#'
#' @keywords internal
#'
#' @importFrom bumphunter clusterMaker
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom IRanges ranges
cluster_sites <- function(methytmle, window_size = 1000) {
  gr <- SummarizedExperiment::rowRanges(methytmle)
  pos <- BiocGenerics::start(IRanges::ranges(gr))
  clusters <- bumphunter::clusterMaker(
    chr = GenomeInfoDb::seqnames(gr),
    pos = pos,
    assumeSorted = FALSE,
    maxGap = window_size
  )
  methytmle@clusters <- as.numeric(clusters)
  return(methytmle)
}

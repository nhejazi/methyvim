library(GenomicRanges)

# build GRanges structure of CpG sites from annotation data
cpg_gr <- makeGRangesFromDataFrame(annot, ignore.strand = TRUE,
                                   seqnames.field = "CHR",
                                   start.field = "MAPINFO",
                                   end.field = "MAPINFO")

# extract annotation data from grouped object and make GRanges
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

# set up and genomic CpG matrix with M-value measurements
clusters <- cluster_sites(cpg_gr)

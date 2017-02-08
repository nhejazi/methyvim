#' constructor for class methadapt (subclass of SummarizedExperiment)
#'
#' @return ".methadapt" S4 object (subclass of "SummarizedExperiment") for DNA
#'         methylation analysis using targeted minimum loss-based estimation.
#'
.methadapt <- setClass(
  "methadapt",
  contains = "SummarizedExperiment",
  slots = list(call = "call")
)

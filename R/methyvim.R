#' Constructor for class methyvim
#'
#' @return class \code{methyvim} object, sub-classed from SummarizedExperiment.
#'
#' @importFrom methods setClass
#' @importFrom tibble data_frame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export .methyvim
#' @exportClass methyvim
#'
.methyvim <- methods::setClass(
       Class = "methyvim",
       slots = list(call = "call",
                    tmleOut = "data_frame",
                    topTable = "data_frame"),
       contains = "SummarizedExperiment"
)

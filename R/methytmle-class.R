#' Constructor for class methytmle
#'
#' @return \code{methytmle} object, subclassed from \code{GenomicRatioSet}.
#'
#' @importFrom methods setClass
#' @importClassesFrom minfi GenomicRatioSet
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'  RangedSummarizedExperiment
#' @import S4Vectors
#'
#' @export .methytmle
#' @exportClass methytmle
#'
#' @examples
#' library(methyvimData)
#' suppressMessages(library(SummarizedExperiment))
#' data(grsExample)
#' # cast the GenomicRatioSet to class methytmle
#' methy_tmle <- .methytmle(grsExample)
#
.methytmle <- methods::setClass(
  Class = "methytmle",
  slots = list(
    call = "call",
    screen_ind = "numeric",
    clusters = "numeric",
    var_int = "numeric",
    param = "character",
    vim = "data.frame",
    ic = "matrix"
  ),
  contains = "GenomicRatioSet"
)

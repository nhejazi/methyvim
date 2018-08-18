#' Constructor for class methytmle
#'
#' @return \code{methytmle} object, subclassed from \code{GenomicRatioSet}.
#'
#' @importFrom methods setClass
#' @importClassesFrom minfi GenomicRatioSet
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'  RangedSummarizedExperiment
#' @importClassesFrom S4Vectors Vector Annotated
#' @import BiocGenerics
#'
#' @export .methytmle
#' @exportClass methytmle
#'
#' @rdname methytmle-class
#'
#' @examples
#' library(methyvimData)
#' suppressMessages(library(SummarizedExperiment))
#' data(grsExample)
#' # cast the GenomicRatioSet to class methytmle
#' methy_tmle <- .methytmle(grsExample)
#' methy_tmle
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

################################################################################

#' @importFrom methods setMethod callNextMethod
#
methods::setMethod("show", "methytmle", function(object) {
  methods::callNextMethod()
  cat("Target Parameter: ")
  cat(param(object))
  cat("\nResults: \n")
  show(vim(object))
})

################################################################################

#' Accessor for Parameter Information
#'
#' @param object S4 object of class \code{methytmle}.
#'
#' @rdname methytmle-class
#'
#' @keywords internal
#
param <- function(object) {
  stopifnot(class(object) == "methytmle")
  object@param
}

################################################################################

#' Accessor for Variable Importance Measure Information
#'
#' @param object S4 object of class \code{methytmle}.
#'
#' @rdname methytmle-class
#'
#' @keywords internal
#'
#' @export
#
vim <- function(object) {
  stopifnot(class(object) == "methytmle")
  object@vim
}


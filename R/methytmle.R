#' Constructor for class methytmle
#'
#' @return \code{methytmle} object, subclassed from \code{GenomicRatioSet}.
#'
#' @importFrom methods setClass
#' @importClassesFrom minfi GenomicRatioSet
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

#' @importFrom methods setGeneric setMethod
#
methods::setGeneric("param", function(object) standardGeneric("param"))
methods::setMethod("param", "methytmle", function(object) {
  object@param
})

################################################################################

#' @importFrom methods setGeneric setMethod
#
methods::setGeneric("vim", function(object) standardGeneric("vim"))
methods::setMethod("vim", "methytmle", function(object) {
  object@vim
})


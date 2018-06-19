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

#' @importFrom methods setMethod
#
methods::setMethod("show", "methytmle", function(object) {
  print(noquote(paste("class:", class(object)[1])))
  print(noquote(paste("data dimension:", paste(as.character(dim(object))[1],
                      as.character(dim(object))[2]))))
  print(noquote(paste("annotation:", object@annotation)))
  print(noquote(paste("target parameter:", object@param)))
  print(noquote("results:"))
  print(object@vim)
})


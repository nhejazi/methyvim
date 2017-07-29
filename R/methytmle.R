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
.methytmle <- methods::setClass(
       Class = "methytmle",
       slots = list(call = "call",
                    g = "matrix",
                    Q = "matrix",
                    ic = "data.frame",
                    vim = "data.frame"),
       contains = "GenomicRatioSet"
)

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
                    screen_ind = "numeric",
                    clusters = "numeric",
                    param = "character",
                    vim = "data.frame",
                    ic = "matrix"),
       contains = "GenomicRatioSet"
)

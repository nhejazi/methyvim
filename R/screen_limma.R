#' Screening procedure based on LIMMA
#'
#' Screens the \code{methytmle} object by way of a \code{limma} model with
#' empirical Bayes shrinkage for CpG sites below the pre-defined p-value cutoff.
#' pre-defined p-value cutoff.
#'
#' ## or something like this? ##
#'
#' Reduces the \code{methytmle} object by way of a \code{limma} model with
#' empirical Bayes shrinkage to include only the CpG sites below the
#' pre-defined p-value cutoff.
#'
#' @param methytmle An object of class \code{methytmle}, with phenotype-level
          matrix missing data dropped if so
#' @param var_int Numeric indicating the column index of the variable of
#'        interest, whether exposure or outcome. If argument \code{vim} is set
#'        to the ATE, then the variable of interest is treated as an exposure;
#'        it is treated as an outcome if this is set to be the NPVI.
#' @param type Character indicating the particular measure of DNA methylation to
#'        be used as the observed data in the estimation procedure, either Beta
#'        values or M-values. The data are accessed via \code{minfi::getBeta} or
#'        \code{minfi::getM}.
#' @param cutoff Numeric indicating the p-value cutoff that defines which
#'        sites pass through the \code{filter}.
#'
#' @return An object of class \code{methytmle}, that is a subset of the original
#' methytmle input object which has been reduced include only the relevant CpG
#' sites. This subset is then used to estimate VIMs.
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom minfi getBeta getM
#' @importFrom SummarizedExperiment colData
#'
limma_screen <- function(methytmle, var_int, type, cutoff = 0.05) {
  # setup design matrix
  design <- as.numeric(SummarizedExperiment::colData(methytmle)[, var_int])
  design <- as.matrix(cbind(rep(1, times = length(design)), design))

  # create expression object for modeling
  if (type == "Beta") {
    methytmle_exprs <- minfi::getBeta(methytmle)
  } else if (type == "Mval") {
    methytmle_exprs <- minfi::getM(methytmle)
  }

  # fit limma model and apply empirical Bayes shrinkage
  mod_fit <- limma::lmFit(object = methytmle_exprs, design = design)
  mod_fit <- limma::eBayes(mod_fit)

  # extract indices of relevant CpG sites
  tt_fit <- limma::topTable(mod_fit, coef = 2, num = Inf, sort.by = "none")
  indices_pass <- which(tt_fit$P.Value < cutoff)

  # add to appropriate slot in the methytmle input object
  methytmle@screen_ind <- as.numeric(indices_pass)
  return(methytmle)
}

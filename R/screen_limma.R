#' Screening procedure based on LIMMA
#'
#' Reduces the \code{methytmle} object by way of a \code{limma} model with
#' empirical Bayes shrinkage to include only the CpG sites below the preset
#' p-value cutoff. INTERNAL USE ONLY.
#'
#' @param methytmle An object of class \code{methytmle}.
#' @param var_int A \code{numeric} vector containing subject-level measurements
#'        of the variable of interest. The length of this vector must match the
#'        number of subjects exactly. If argument \code{vim} is set to "ate" or
#'        "rr", then the variable of interest is treated as an exposure, and the
#'        variable must be binary in such cases. If setting \code{vim} to target
#'        parameters assessing continuous treatment effects, then the variable
#'        need not be binary of course.
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
#' @keywords internal
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom minfi getBeta getM
#
limma_screen <- function(methytmle, var_int, type, cutoff = 0.05) {

  stopifnot(class(methytmle) == "methytmle")

  # setup design matrix
  design <- as.matrix(cbind(rep(1, times = length(var_int)), var_int))

  # create expression object for modeling
  # M-values have the property of having far more centrality than Beta-values.
  # Parametric techniques work best with an underlying Gaussian distribution to
  # the data so M values might be more appropriate to use for analysis.
  if (type == "Beta") {
    Betas <- minfi::getBeta(methytmle)
    # shifting the betas from 0 and 1 prevents Inf/Nan Mvals
    newBetas <- Harman::shiftBetas(methytmle_exprs, shiftBy=1e-4)
    methytmle_exprs <- wateRmelon::Beta2M(newBetas)
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

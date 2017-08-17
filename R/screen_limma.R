#' Screening procedure based on LIMMA
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

#' Screening procedure based on LIMMA
#'
#' @importFrom limma lmFit eBayes topTable
#'
#' @export lmfit_screen
#'

lmfit_screen <- function(meth, design, cutoff, ...) {
  mod_fit <- limma::lmFit(meth, design)
  mod_fit <- limma::eBayes(mod_fit)
  tt_out <- limma::topTable(mod_fit, coef = 2, num = Inf, sort.by = "none")
  indices_pass <- which(tt_out$P.Value < cutoff)
  tt_screened <- subset(tt_out, P.Value < cutoff)
  out <- list(indices_pass, tt_screened)
  return(out)
}

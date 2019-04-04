#' FDR-MSA correction
#'
#' Modified FDR Controlling Procedure for Multi-Stage Analyses (MJ van der Laan
#' and C Tuglus, 2009, <doi:10.2202/1544-6115.1397>)
#'
#' @param pvals Numeric vector containing the p-values that result from any
#'  chosen statistical hypothesis testing procedure.
#' @param total_obs Numeric indicating the total number of observations that
#'  would have been available for testing prior to the selection procedure
#'  employed in the multi-stage analysis performed.
#'
#' @return A \code{numeric} vector of corrected p-values, controlling the False
#'  Discovery Rate, using the method of Tuglus and van der Laan.
#'
#' @importFrom stats p.adjust
#'
#' @export
#'
#' @examples
#' g <- 1e4
#' n <- 1e2
#' p <- abs(rnorm(n, mean = 1e-8, sd = 1e-2))
#' # treating the vector p as one of p-values, FDR-MSA may be applied
#' fdr_p <- fdr_msa(pvals = p, total_obs = g)
fdr_msa <- function(pvals, total_obs) {
  pvals_not_tested <- rep(1, total_obs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fdr_adj <- stats::p.adjust(pvals_all, method = "fdr")
  fdr_out <- fdr_adj[seq_along(pvals)]
  return(fdr_out)
}

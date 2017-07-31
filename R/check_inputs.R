#' Check Function Inputs
#'
#' Check the input values of function parameters for errors.
#'
#' ...
#' ...
#'
#' @return Options to be passed to various \code{tmle_*} functions.
#'

check_inputs <- function(catch_inputs) {
  message("not yet implemented")
  # NOTE: make sure to check combinations of the variable of interest and the
  # stated parameter of interest (i.e., VIM). In particular, we will use the
  # following heuristic to decide between NPVI and ATE: if the variable of
  # interest is discretized, use the ATE, with the Y = CpG site i, A = treatment
  # (discretized), and W = neighboring uncorrelated sites; if the variable of
  # interest is continuous, use NPVI. We don't actually need to impose the
  # heuristic, but use it to check that the user has specified the ATE vs. NPVI
  # parameters correctly.
}

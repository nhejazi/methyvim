#' Check Function Inputs
#'
#' Check the input values of function parameters for errors.
#'
#' @param data ...
#' @param var ...
#' @param vim ...
#' @param type ..
#' @param filter ...
#' @param filter_cutoff ...
#' @param window ...
#' @param corr ...
#' @param obs_per_covar ...
#' @param parallel ...
#' @param future ...
#' @param bppar ...
#' @param return_ic ...
#' @param shrink_ic ..
#' @param tmle_type ...
#' @param tmle_args ...
#'
#' @return Options to be passed to various \code{tmle_*} functions.
#'
check_inputs <- function(data, var, vim, type,
                         filter, filter_cutoff, window, corr, obs_per_covar,
                         parallel, future, bppar,
                         return_ic, shrink_ic,
                         tmle_type, tmle_args) {

  message("This function is only partially implemented currently.")
  # NOTE: make sure to check combinations of the variable of interest and the
  # stated parameter of interest (i.e., VIM). In particular, we will use the
  # following heuristic to decide between NPVI and ATE: if the variable of
  # interest is discretized, use the ATE, with the Y = CpG site i, A = treatment
  # (discretized), and W = neighboring uncorrelated sites; if the variable of
  # interest is continuous, use NPVI. We don't actually need to impose the
  # heuristic, but use it to check that the user has specified the ATE vs. NPVI
  # parameters correctly.

  if (vim == "npvi" & is.null(tmle_args$npvi_descr)) {
    npvi_descr_defaults <- list(f = identity, iter = 10, cvControl = 2,
                                nMax = 30,
                                stoppingCriteria = list(mic = 0.001,
                                                        div = 0.001,
                                                        psi = 0.01)
                               )
    tmle_args$npvi_descr <- npvi_descr_defaults
  }

  if (tmle_type == "glm") {
    if (vim == "ate") {
      # set GLM libraries for "tmle" package
      g_lib <- "SL.mean"
      Q_lib <- "SL.glm"
      tmle_args$g_lib <- g_lib
      tmle_args$Q_lib <- Q_lib
    } else if (vim == "npvi") {
      message("finding a way to select GLM flavor for tmle.npvi")
    }
  }

  #  TODO: all other arguments in catch_inputs
}

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
  message("check_inputs is only partially implemented currently.")
  # NOTE: make sure to check combinations of the variable of interest and the
  # stated parameter of interest (i.e., VIM). In particular, we will use the
  # following heuristic to decide between NPVI and ATE: if the variable of
  # interest is discretized, use the ATE, with the Y = CpG site i, A = treatment
  # (discretized), and W = neighboring uncorrelated sites; if the variable of
  # interest is continuous, use NPVI. We don't actually need to impose the
  # heuristic, but use it to check that the user has specified the ATE vs. NPVI
  # parameters correctly.

  if (catch_inputs$vim == "NPVI" & is.null(catch_inputs$tmle_args$npvi_descr)) {
    npvi_descr_defaults <- list(f = identity, iter = 10, cvControl = 2,
                                nMax = 30,
                                stoppingCriteria = list(mic = 0.001,
                                                        div = 0.001,
                                                        psi = 0.01)
                               )
    catch_inputs$tmle_args$npvi_descr <- npvi_descr_defaults
  }

  if (catch_inputs$tmle_type == "glm") {
    if (catch_inputs$vim == "ATE") {
      # set GLM libraries for "tmle" package
      g_lib <- "SL.mean"
      Q_lib <- "SL.glm"
      catch_inputs$tmle_args$g_lib <- g_lib
      catch_inputs$tmle_args$Q_lib <- Q_lib
    } else if (catch_inputs$vim == "NPVI") {
      message("finding a way to select GLM flavor for tmle.npvi")
    }
  }

  #  TODO: all other arguments in catch_inputs
}

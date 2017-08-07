#' Differential Methylation Statistics with Variable Importance Measures
#'
#' Computes the Targeted Minimum Loss-Based Estimate of a specified statistical
#' target parameter, formally defined within models from causal inference. The
#' variable importance measures currently supported are the Average Treatment
#' Effect (ATE) and a Nonparametric Variable Importance Measure (NPVI, formally
#' defined by Chambaz, Neuvial, and van der Laan <doi:10.1214/12-EJS703>).
#'
#' @param data_grs An object of class \code{minfi::GenomicRatioSet}, containing
#'        standard data structures associated with DNA Methylation experiments.
#'        Consult the documentation for \code{minfi} to construct such objects.
#' @param var_int Numeric indicating the column index of the variable of
#'        interest, whether exposure or outcome. If argument \code{vim} is set
#'        to the ATE, then the variable of interest is treated as an exposure;
#'        it is treated as an outcome if this is set to be the NPVI.
#' @param vim Character indicating the variable importance measure to be used
#'        in the estimation procedure. Currently supported options are the ATE
#'        for discretized exposures and NPVI for continuous exposures. ATE is
#'        the appropriate choice when the underlying scientific question is of
#'        the effect of an exposure on methylation, while NPVI is the parameter
#'        of choice when the effect of methylation on an outcome is sought.
#' @param type Character indicating the particular measure of DNA methylation to
#'        be used as the observed data in the estimation procedure, either Beta
#'        values or M-values. The data are accessed via \code{minfi::getBeta} or
#'        \code{minfi::getM}.
#' @param filter Character indicating...
#'        ...
#' @param filter_cutoff Numeric indicating...
#'        ...
#' @param window_bp Numeric indicating...
#'        ...
#' @param corr_max Numeric indicating...
#'        ...
#' @param obs_per_covar Numeric indicating...
#'        ...
#' @param parallel Logical or Numeric indicating...
#'        ...
#' @param return_ic Logical indicating...
#'        ...
#' @param shrink_ic Logical indicating...
#'        ...
#' @param tmle_type Character indicating the general class of regression models
#'        to be used in fitting the propensity score and outcome regressions.
#'        This is generally a shorthand and is overridden by \code{tmle_args} if
#'        that argument is changed from its default values.
#' @param tmle_args List giving several key arguments to be passed to one of
#'        \code{tmle::tmle} or \code{tmle.npvi::tmle.npvi}, depending on the
#'        particular variable importance measure specified. This overrides
#'        \code{tmle_type}, which itself provides sensible defaults. Consider
#'        changing this away from default settings only if you have sufficient
#'        experience with theory and software for targeted learning. For more
#'        information, consider consulting the documentation of the \code{tmle}
#'        and \code{tmle.npvi} packages.
#'
#' @return An object of class \code{methytmle}, with all unique slots filled in,
#'         in particular, including indices of CpG sites that pass screening,
#'         cluster of neighboring CpG sites, and a matrix of the results of the
#'         estimation procedure performed for the given variable importance
#'         measure. Optionally, estimates of the propensity score and outcome
#'         regressions, as well as the original data rotated into influence
#'         curve space may be returned, if so requested.
#'
#' @export methyvim
#'

methyvim <- function(data_grs,
                     var_int = 1,
                     vim = c("ATE", "NPVI"),
                     type = c("Beta", "Mval"),
                     filter = c("limma", "npvi", "adaptest"),
                     filter_cutoff = 0.05,
                     window_bp = 1e3,
                     corr_max = 0.35,
                     obs_per_covar = 15,
                     parallel = TRUE,
                     return_ic = FALSE,
                     shrink_ic = FALSE,
                     tmle_type = c("glm", "super_learning"),
                     tmle_args = list(family = "binomial",
                                      g_lib = c("SL.mean", "SL.glm",
                                                "SL.randomForest"),
                                      Q_lib = c("SL.mean", "SL.glm"),
                                      npvi_cutoff = 0.25,
                                      npvi_descr = NULL)
                    ) {

  # ============================================================================
  # catch input and return in output object for user convenience
  # ============================================================================
  call <- match.call(expand.dots = TRUE)
  filter <- match.arg(filter)
  type <- match.arg(type)
  vim <- match.arg(vim)

  # ============================================================================
  # catch inputs to pass to downstream functions (for estimation and such)
  # ============================================================================
  catch_inputs <- list(data = data_grs,
                       var = var_int,
                       vim = vim,
                       type = type,
                       filter = filter,
                       filter_cutoff = filter_cutoff,
                       window = window_bp,
                       corr = corr_max,
                       obs_per_covar = obs_per_covar,
                       par = parallel,
                       return_ic = return_ic,
                       shrink_ic = shrink_ic,
                       tmle_type = tmle_type,
                       tmle_args = tmle_args)
  # check that inputs satisfy expectations
  ## and clean up check_inputs (e.g., rm NPVI stuff if ATE is specified)
  catch_inputs <- check_inputs(catch_inputs)

  # ============================================================================
  # invoke S4 class constructor for "methadapt" object
  # ============================================================================
  methy_tmle <- .methytmle(catch_inputs$data)
  methy_tmle@call <- call

  #=============================================================================
  # set up parallelization if so desired
  # ============================================================================
  if (catch_inputs$par) {
    set_parallel()
  }

  # ============================================================================
  # operate on the type of data specified
  # ============================================================================
  #if (catch_inputs$type == "Beta") {
  #  extract_measures <- parse(text = "getBeta(methy_tmle)")
  #} else if (catch_inputs$type == "Mval") {
  #  extract_measures <- parse(text = "getM(methy_tmle)")
  #}

  #=============================================================================
  # check if there is missing data in the phenotype-level matrix and drop if so
  # ============================================================================
  cols_na <- colSums(sapply(colData(methy_tmle), is.na))
  na_var_id <- names(which(cols_na != 0))
  if (length(na_var_id) != 0) {
    message(paste("Missing data detected: Dropping variable(s)", na_var_id,
                  "from phenotype-level data matrix."))
  }
  colData(methy_tmle)[na_var_id] <- NULL

  #=============================================================================
  # screen sites to produce a subset on which to estimate VIMs
  # ============================================================================
  if (catch_inputs$filter == "limma") {
    methy_tmle <- limma_screen(methytmle = methy_tmle,
                               var_int = catch_inputs$var,
                               type = catch_inputs$type)
  } else {
    stop("The designated filtering method has not yet been implemented.")
  }

  #=============================================================================
  # cluster sites based on genomic windows
  #=============================================================================
  methy_tmle <- cluster_sites(methy_tmle = methy_tmle)

  #=============================================================================
  # TMLE procedure for targeted differential methylation analysis
  # ============================================================================
  if (catch_inputs$vim == "ATE") {
    methy_tmle@vim <- methyvim_ate(methy_tmle = methy_tmle,
                                   catch_inputs_ate = catch_inputs)
  } else if (catch_inputs$vim == "NPVI") {
    methy_tmle@vim <- methyvim_npvi(methy_tmle = methy_tmle,
                                    catch_inputs_npvi = catch_inputs)
  } else {
    stop("The specified variable importance parameter is not available.")
  }

  # NOTE: what else do we do before returning output...
}

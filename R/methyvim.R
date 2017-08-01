#' Differential Methylation Analysis with Nonparametric Variable Importance
#'
#' Computes the Targeted Minimum Loss-Based Estimate of...
#'
#' @param data_grs - ...
#'        ...
#' @param var_int - ...
#'        ...
#' @param min_sites Numeric indicating the minimum number of sites to be kept
#'        after the application of a filtering procedure.
#' @param neighbors Numeric indicating the size (in bp) of genomic windows to
#'        use when assigning CpG sites to groups.
#' @param parallel c(logical, numeric) - whether to use or number of cores to be
#'        used in parallelizing the TMLE-based procedure across genomic sites.
#' @param family (character) - specification of error family; available options
#'        are "binomial" or "gaussian" (default).
#' @param g_lib - library of learning algorithms to be used in fitting the model
#'        to "g" in the Targeted Minimum Loss-Based Estimation procedure.
#' @param Q_lib - library of learning algorithms to be used in fitting the model
#'        to "Q" in the Targeted Minimum Loss-Based Estimation procedure.
#'
#' @export methyvim
#'

methyvim <- function(data_grs,
                     var_int = 1,
                     cpg_is = "exposure",
                     type = c("Beta", "Mval"),
                     vim = c("ATE", "NPVI"),
                     filter = c("limma", "npvi", "adaptest"),
                     neighbors = 1e3,
                     corr_max = 0.35,
                     preprocess = NULL,
                     parallel = TRUE,
                     dimen_red = FALSE,
                     return_ic = FALSE,
                     shrink_ic = FALSE,
                     family = "gaussian",
                     g_lib = c("SL.mean", "SL.glm", "SL.randomForest"),
                     Q_lib = c("SL.mean", "SL.randomForest"),
                     subj_per_covar = 15
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
  catch_inputs <- list(data = data_grs, var = var_int, cpg_iss = cpg_is,
                       type = type, vim = vim, neighbors = neighbors,
                       preprocess = preprocess, filter = filter,
                       min_sites = min_sites, family = family,
                       g_lib = g_lib, Q_lib = Q_lib, parallel = parallel,
                       return_ic = return_ic, shrink_ic = shrink_ic,
                       dm = dimen_red, corr_max = corr_max,
                       subj_per_covar = subj_per_covar)
  catch_inputs <- check_inputs(catch_inputs)

  # ============================================================================
  # invoke S4 class constructor for "methadapt" object
  # ============================================================================
  methy_tmle <- .methytmle(catch_inputs$data)
  methy_tmle@call <- call

  #=============================================================================
  # set up parallelization if so desired
  # ============================================================================
  if (catch_inputs$parallel) {
    set_parallel()
  }

  # ============================================================================
  # operate on the type of data specified
  # ============================================================================
  if (catch_inputs$type == "Beta") {
    extract_measures <- parse(text = "getBeta(methy_tmle)")
  } else if (catch_inputs$type == "Mval") {
    extract_measures <- parse(text = "getM(methy_tmle)")
  }

  # ============================================================================
  # preprocess array data if so requested
  # ============================================================================
  if (!is.null(catch_inputs$preprocess)) {
    # methy_tmle <- do_preprocess(methy_tmle)
    message("Support for preprocesing is planned but not yet implemented.")
  }

  #=============================================================================
  # check if there is missing data in the phenotype-level matrix and drop if so
  # ============================================================================
  cols_na <- colSums(sapply(colData(methy_tmle), is.na))
  na_var_id <- names(which(cols_na != 0))
  if (length(na_var_id) != 0) {
    message(paste("Missing data detected: Dropping variable(s)", na_var_id,
                  "from phenotype matrix."))
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
    stop("The chosen filtering method has not yet been implemented.")
  }

  #=============================================================================
  # cluster sites based on genomic windows
  #=============================================================================
  methy_tmle <- cluster_sites(methy_tmle = methy_tmle)

  #=============================================================================
  # TMLE procedure for targeted differential methylation analysis
  # ============================================================================
  if (catch_inputs$vim == "ATE") {
    methy_tmle <- methyvim_ate(methy_tmle)
  } else if (catch_inputs$vim == "NPVI") {
    methy_tmle <- methyvim_npvi(methy_tmle)
  } else {
    stop("The specified variable importance parameter is not available.")
  }

  # NOTE: what else do we do before returning output...
}

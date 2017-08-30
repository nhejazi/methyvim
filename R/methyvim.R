utils::globalVariables(c("colData<-"))

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
#' @param sites_comp Numeric indicating the number of sites for which a variable
#'        importance measure is to be estimated post-screening.
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
#' @param filter Character indicating the model to be implemented when screening
#'        the \code{data_grs} object for CpG sites. Currently
#'        supported options are limma, npvi, and adaptest.  (?)
#' @param filter_cutoff Numeric indicating the p-value cutoff that defines which
#'        sites pass through the \code{filter}.
#' @param window_bp Numeric indicating the maximum distance between two sites
#'        for them to be considered neighboring sites (?)
#' @param corr_max Numeric indicating the maximum correlation that a neighboring
#'        site can have with the target site.
#' @param obs_per_covar Numeric indicating the number of observations needed for
#'        for covariate included in W for downstream analysis. This ensures the
#'        data is sufficient to control for the covariates.
#' @param parallel Logical indicating whether parallelization ought to be used.
#'        See the documentation of \code{set_parallel} for more information, as
#'        this arugment is passed directly to that internal function.
#' @param future_param Character indicating the type of parallelization to be
#'        used from the list available via the \code{future} package. See the
#'        documentation for \code{set_parallel} for more information, as this
#'        argument is passed directly to that internal function.
#' @param bppar_type Character specifying the type of backend to be used for
#'        parallelization via \code{BiocParallel}. See the documentation for
#'        \code{set_parallel} for more information, as this argument is passed
#'        directly to that internal function.
#' @param return_ic Logical indicating whether an influence curve estimate
#'          should be returned for each site that passed through the filter.
#' @param shrink_ic Logical indicating whether limma should be applied to reduce
#'        the variance in the ic based estimates in \code{return_ic}.
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
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel register bpprogressbar DoparParam
#' @importFrom future plan multiprocess sequential
#' @importFrom doFuture registerDoFuture
#'
#' @export
#'
methyvim <- function(data_grs,
                     sites_comp = 10,
                     var_int = 1,
                     vim = c("ate", "rr", "npvi"),
                     type = c("Beta", "Mval"),
                     filter = c("limma"),
                     filter_cutoff = 0.05,
                     window_bp = 1e3,
                     corr_max = 0.75,
                     obs_per_covar = 20,
                     parallel = TRUE,
                     future_param = NULL,
                     bppar_type = NULL,
                     return_ic = FALSE,
                     shrink_ic = FALSE,
                     tmle_type = c("glm", "sl"),
                     tmle_args = list(family = "binomial",
                                      g_lib = NULL, Q_lib = NULL,
                                      npvi_cutoff = 0.25, npvi_descr = NULL)
                    ) {
  # ============================================================================
  # catch input for user convenience and check input types where possible
  # ============================================================================
  # catch function call
  call <- match.call(expand.dots = TRUE)
  # type checking (in same order the arguments appear)
  vim <- match.arg(vim)
  type <- match.arg(type)
  filter <- match.arg(filter)
  tmle_type <- match.arg(tmle_type)

  # ============================================================================
  # set treatment mechanism and outcome regression libraries
  # ============================================================================
  if(is.null(tmle_args$g_lib)) {
    tmle_args$g_lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  }
  if(is.null(tmle_args$Q_lib)) {
    tmle_args$Q_lib <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  }

  # ============================================================================
  # modify arguments to TMLE functions based on type of type requested
  # ============================================================================
  if (tmle_type == "glm") {
    if (vim %in% c("ate", "rr")) {
      # set GLM libraries for "tmle" package
      tmle_args$g_lib <- c("SL.mean", "SL.glm")
      tmle_args$Q_lib <- "SL.glm"
    }
  }

  # ============================================================================
  # if NPVI parameter requested, set sensible defaults
  # ============================================================================
  if (vim == "npvi" & is.null(tmle_args$npvi_descr)) {
    npvi_descr_defaults <- list(f = identity, iter = 10, cvControl = 2,
                                nMax = 30,
                                stoppingCriteria = list(mic = 0.001,
                                                        div = 0.001,
                                                        psi = 0.01)
                               )
    tmle_args$npvi_descr <- npvi_descr_defaults
  }

  # ============================================================================
  # invoke S4 class constructor for object
  # ============================================================================
  methy_tmle <- .methytmle(data_grs)
  methy_tmle@call <- call
  methy_tmle@var_int <- var_int

  #=============================================================================
  # set up parallelization if so desired
  # ============================================================================
  set_parallel(parallel = parallel,
               future_param = future_param,
               bppar_type = bppar_type)

  #=============================================================================
  # check if there is missing data in the phenotype-level matrix and drop if so
  # ============================================================================
  cols_na <- colSums(sapply(SummarizedExperiment::colData(methy_tmle), is.na))
  na_var_id <- names(which(cols_na != 0))
  if (length(na_var_id) != 0) {
    message(paste("Missing data detected: Dropping variable(s)", na_var_id,
                  "from phenotype-level data matrix."))
  }
  colData(methy_tmle)[na_var_id] <- NULL

  #=============================================================================
  # screen sites to produce a subset on which to estimate VIMs
  # ============================================================================
  if (filter == "limma") {
    methy_tmle <- limma_screen(methytmle = methy_tmle,
                               var_int = var_int,
                               type = type)
  }

  #=============================================================================
  # cluster sites based on genomic windows
  #=============================================================================
  methy_tmle <- cluster_sites(methytmle = methy_tmle)

  #=============================================================================
  # TMLEs for the Average Treatment Effect (ATE) and Risk Ratio (RR) parameters
  # ============================================================================
  if (vim %in% c("ate", "rr")) {
    # make sure that the outcome data is of class numeric
    var_of_interest <- SummarizedExperiment::colData(methy_tmle)[, var_int]
    var_of_interest <- as.numeric(var_of_interest)
    if (class(var_of_interest) != "numeric") {
      var_of_interest <- as.numeric(var_of_interest)
    }

    # get names of sites to be added to output object
    cpg_screened_names <- names(methy_tmle[methy_tmle@screen_ind])

    ## TODO: THIS IS FOR TESTING ONLY
    cpg_screened_names <- cpg_screened_names[seq_len(sites_comp)]

    # object of screened CpG site indices to loop over in TMLE procedure
    methy_tmle_ind <- seq_along(methy_tmle@screen_ind)

    methy_vim_out <- BiocParallel::bplapply(methy_tmle_ind[seq_len(sites_comp)],
                                            FUN = methyvim_tmle,
                                            methytmle_screened = methy_tmle,
                                            var_of_interest = var_of_interest,
                                            type = type,
                                            corr = corr_max,
                                            obs_per_covar = obs_per_covar,
                                            target_param = vim,
                                            g_lib = tmle_args$g_lib,
                                            Q_lib = tmle_args$Q_lib,
                                            family = tmle_args$family,
                                            return_ic = return_ic
                                           )
    methy_vim_out <- do.call(rbind.data.frame, methy_vim_out)

    # TMLE procedure is now done, so let's just make the output object pretty...
    if (vim == "ate") {
      colnames(methy_vim_out) <- c("lowerCI_ATE", "est_ATE", "upperCI_ATE",
                                   "Var_ATE", "pval", "n_neighbors_all",
                                   "n_neighbors_w", "max_corr_w")
      methy_tmle@param <- "Average Treatment Effect"
    } else {
      colnames(methy_vim_out) <- c("lowerCI_logRR", "est_logRR",
                                   "upperCI_logRR", "Var_logRR", "pval",
                                   "n_neighbors_all", "n_neighbors_w",
                                   "max_corr_w")
      methy_tmle@param <- "Risk Ratio"
    }

    rownames(methy_vim_out) <- cpg_screened_names
    methy_tmle@vim <- methy_vim_out
  }
  # Let's give 'em some output
  return(methy_tmle)
}

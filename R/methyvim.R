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
#'        for covariate included in (W/downstream analysis). This ensures the
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
#' @importFrom BiocParallel bplapply
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
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
                     future_param = NULL,
                     bppar_type = NULL,
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
# if you say glm it will choose mean and glm
# if you choose super learner it will choose more algorithms
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
  # check that inputs satisfy expectations
  ## and clean up check_inputs (e.g., rm NPVI stuff if ATE is specified)
  #check_inputs(data = data_grs,
  #             var = var_int,
  #             vim = vim,
  #             type = type,
  #             filter = filter,
  #             filter_cutoff = filter_cutoff,
  #             window = window_bp,
  #             corr = corr_max,
  #             obs_per_covar = obs_per_covar,
  #             parallel = parallel,
  #             future = future_param,
  #             bppar = bppar_type,
  #             return_ic = return_ic,
  #             shrink_ic = shrink_ic,
  #             tmle_type = tmle_type,
  #             tmle_args = tmle_args
  #            )

  # ============================================================================
  # invoke S4 class constructor for "methadapt" object
  # ============================================================================
  methy_tmle <- .methytmle(data)
  methy_tmle@call <- call

  #=============================================================================
  # set up parallelization if so desired
  # ============================================================================
  #set_parallel(parallel = parallel,
  #             future_param = future_param,
  #             bppar_type = bppar_type)
  n_cores <- parallel::detectCores()
  doParallel::registerDoParallel(n_cores)

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
  if (filter == "limma") {
    methy_tmle <- limma_screen(methytmle = methy_tmle,
                               var_int = var_int,
                               type = type)
  } else {
    stop("The designated filtering method has not yet been implemented.")
  }

  #=============================================================================
  # cluster sites based on genomic windows
  #=============================================================================
  methy_tmle <- cluster_sites(methytmle = methy_tmle)

  #=============================================================================
  # ATE TMLE procedure for targeted differential methylation analysis
  # ============================================================================
  if (vim == "ATE") {

    # make sure that the outcome data is of class numeric
    var_of_interest <- as.numeric(colData(methy_tmle)[, var_int])
    if (class(var_of_interest) != "numeric") {
      var_of_interest <- as.numeric(var_of_interest)
    }

    # get names of sites to be added to output object
    cpg_screened_names <- names(methy_tmle[methy_tmle@screen_ind])
    cpg_screened_names <- cpg_screened_names[seq_len(100)] ## TODO: REMOVE, FOR TESTING ONLY

    # object of screened CpG site indices to loop over in TMLE procedure
    methy_tmle_ind <- seq_along(methy_tmle@screen_ind)
    methy_tmle_ind <- methy_tmle_ind[seq_len(100)] ## TODO: REMOVE, FOR TESTING
    #methy_vim_out <- BiocParallel::bplapply(X = methy_tmle_ind,
    #                                        FUN = methyvim_ate,
    #                                        methy_tmle_screened = methy_tmle,
    #                                        var_of_interest = var_of_interest,
    #                                        type = "Mval",
    #                                        corr = 0.80,
    #                                        obs_per_covar = 20,
    #                                        family = "gaussian",
    #                                        g_lib = c("SL.mean", "SL.glm"),
    #                                        Q_lib = c("SL.mean", "SL.glm"),
    #                                        return_ic = FALSE
    #                                       )
    methy_vim_out <- foreach::foreach(i_site = methy_tmle_ind,
                                      .packages = c("tmle"),
                                      .combine = rbind) %dopar% {

      message(paste("Computing targeted estimate for site", i_site, "of",
                    length(methy_tmle_ind)))

      out <- methyvim_ate(target_site = i_site,
                          methy_tmle_screened = methy_tmle,
                          var_of_interest = var_of_interest,
                          type = "Mval",
                          corr = 0.80,
                          obs_per_covar = 20,
                          g_lib = c("SL.mean", "SL.glm"),
                          Q_lib = c("SL.mean", "SL.glm"),
                          family = "gaussian",
                          return_ic = FALSE
                         )
    }
    methy_vim_out <- as.data.frame(methy_vim_out)

    # TMLE procedure is now done, so let's just make the output object pretty...
    #methy_vim_out <- do.call(rbind.data.frame, methy_vim_out)
    colnames(methy_vim_out) <- c("lower_CI_ATE", "est_ATE", "upper_CI_ATE",
                                 "Var", "pval", "n_neighbors_all",
                                 "n_neighbors_w", "max_corr_w")
    rownames(methy_vim_out) <- cpg_screened_names
    methy_tmle@vim <- methy_vim_out

  #=============================================================================
  # NPVI TMLE procedure for targeted differential methylation analysis
  # ============================================================================
  } else if (vim == "NPVI") {
    stop("Support for TMLE-NPVI is planned but not yet implemented.")
    #methy_tmle@vim <- methyvim_npvi(methy_tmle = methy_tmle,
    #                                catch_inputs_npvi = catch_inputs)
  } else {
    stop("The specified variable importance parameter is not available.")
  }

  # NOTE: what else do we do before returning output...
}

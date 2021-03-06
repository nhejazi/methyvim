utils::globalVariables(c("colData<-"))

#' Differential Methylation Statistics with Variable Importance Measures
#'
#' Computes the Targeted Minimum Loss Estimate of a specified statistical target
#' parameter, formally defined within models from causal inference. The variable
#' importance measures currently supported are the Average Treatment Effect
#' (ATE) and a Nonparametric Variable Importance Measure (NPVI, formally
#' defined by Chambaz, Neuvial, and van der Laan <doi:10.1214/12-EJS703>).
#'
#' @param data_grs An object of class \code{\link[minfi]{GenomicRatioSet}},
#'  containing standard data structures for DNA Methylation experiments.
#'  Consult the documentation of \pkg{minfi} to construct such objects.
#' @param var_int A \code{numeric} vector containing subject-level measurements
#'  of the variable of interest. The length of this vector must match the
#'  number of subjects exactly. If argument \code{vim} is set to "ate" or "rr",
#'  then the variable of interest is treated as an exposure, and the variable
#'  must be binary in such cases. If setting \code{vim} to target parameters
#'  assessing continuous treatment effects, then the variable need not be binary
#'  of course.
#' @param vim Character indicating the variable importance measure to be used in
#'  the estimation procedure. Currently supported options are the ATE for
#'  discretized exposures and NPVI for continuous exposures. ATE and RR are the
#'  appropriate choices when the underlying scientific question is of the effect
#'  of an exposure on methylation, while NPVI (and other continuous treatment
#'  parameters) ought to be used when the effect of methylation on an outcome is
#'  sought.
#' @param type Character indicating the particular measure of DNA methylation to
#'  be used as the observed data in the estimation procedure, either Beta values
#'  or M-values. The data are accessed via \code{\link[minfi]{getBeta}} or
#'  \code{\link[minfi]{getM}}.
#' @param filter Character indicating the model to be implemented when screening
#'  the \code{data_grs} object for CpG sites. The only currently supported
#'  option is "limma".
#' @param filter_cutoff Numeric indicating the p-value cutoff that defines which
#'  sites pass through the \code{filter}.
#' @param window_bp Numeric indicating the maximum genomic distance (in base
#'  pairs) between two sites for them to be considered neighboring sites.
#' @param corr_max Numeric indicating the maximum correlation that a neighboring
#'  site can have with the target site.
#' @param obs_per_covar Numeric indicating the number of observations needed for
#'  for covariate included in W for downstream analysis. This ensures the data
#'  is sufficient to control for the covariates.
#' @param sites_comp A \code{numeric} indicating the maximum number of sites for
#'  which a variable importance measure is to be estimated post-screening. This
#'  is not typically useful in scientific settings, but may be useful when a
#'  large number of CpG sites pass the initial screening phase.
#' @param parallel Logical indicating whether parallelization ought to be used.
#'  See the documentation of \code{set_parallel} for more information, as this
#'  argument is passed directly to that internal function.
#' @param future_param Character indicating the type of parallelization to be
#'  used from the list available via the \pkg{future} package. See the
#'  documentation for \code{set_parallel} for more information, as this argument
#'  is passed directly to that internal function.
#' @param bppar_type Character specifying the type of backend to be used for
#'  parallelization via \pkg{BiocParallel}. See the documentation for
#'  \code{set_parallel} for more information, as this argument is passed
#'  directly to that internal function.
#' @param return_ic Logical indicating whether an influence curve estimate
#'  should be returned for each site that passed through the filter.
#' @param shrink_ic Logical indicating whether limma should be applied to reduce
#'  the variance in the ic based estimates in \code{return_ic}.
#' @param tmle_type Character indicating the general class of regression models
#'  to be used in fitting the propensity score and outcome regressions. This is
#'  generally a shorthand and is overridden by \code{tmle_args} if that argument
#'  is changed from its default values.
#' @param tmle_args List giving several key arguments to be passed to one of
#'  \code{\link[tmle]{tmle}} or \code{\link[tmle.npvi]{tmle.npvi}}, depending on
#'  the particular variable importance measure specified. This overrides
#'  \code{tmle_type}, which itself provides sensible defaults. Consider changing
#'  this away from default settings only if you have sufficient experience with
#'  software and the underlying theory for Targeted Learning. For more
#'  information, consider consulting the documentation of the \pkg{tmle} and
#'  \pkg{tmle.npvi} packages.
#' @param tmle_backend A \code{character} indicating the package to be used in
#'  the estimation procedure. The user should only set this parameter if they
#'  have sufficient familiarity with the backend packages used for estimation.
#'  Current choices include \pkg{tmle}, \pkg{drtmle}, and \pkg{tmle.npvi}.
#'
#' @return An object of class \code{methytmle}, with all unique slots filled in,
#'  in particular, including indices of CpG sites that pass screening, cluster
#'  of neighboring CpG sites, and a matrix of the results of the estimation
#'  procedure performed for the given variable importance measure. Optionally,
#'  estimates of the propensity score and outcome regressions, as well as the
#'  original data rotated into influence curve space may be returned, if so
#'  requested.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel register bpprogressbar DoparParam
#' @importFrom future plan multiprocess sequential
#' @importFrom doFuture registerDoFuture
#'
#' @export
#'
#' @examples
#' library(methyvimData)
#' suppressMessages(library(SummarizedExperiment))
#' data(grsExample)
#' var_int <- colData(grsExample)[, 1]
#' # TMLE procedure for the ATE parameter over M-values with Limma filtering
#' methyvim_out_ate <- suppressWarnings(
#'   methyvim(
#'     data_grs = grsExample, sites_comp = 1, var_int = var_int,
#'     vim = "ate", type = "Mval", filter = "limma", filter_cutoff = 0.05,
#'     parallel = FALSE, tmle_type = "sl"
#'   )
#' )
methyvim <- function(data_grs,
                     var_int,
                     vim = c("ate", "rr", "npvi"),
                     type = c("Beta", "Mval"),
                     filter = c("limma"),
                     filter_cutoff = 0.05,
                     window_bp = 1e3,
                     corr_max = 0.75,
                     obs_per_covar = 20,
                     sites_comp = NULL,
                     parallel = TRUE,
                     future_param = NULL,
                     bppar_type = NULL,
                     return_ic = FALSE,
                     shrink_ic = FALSE,
                     tmle_type = c("glm", "sl"),
                     tmle_args = list(
                       g_lib = c("SL.mean", "SL.glm", "SL.bayesglm", "SL.gam"),
                       Q_lib = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                       cv_folds = 5,
                       npvi_cutoff = 0.25,
                       npvi_descr = NULL
                     ),
                     tmle_backend = c("tmle", "drtmle", "tmle.npvi")) {
  # ===========================================================================
  # catch input for user convenience and check input types where possible
  # ===========================================================================
  call <- match.call(expand.dots = TRUE)
  vim <- match.arg(vim)
  type <- match.arg(type)
  filter <- match.arg(filter)
  tmle_type <- match.arg(tmle_type)
  tmle_backend <- match.arg(tmle_backend)

  # check that variable of interest is the correct length
  if (length(var_int) != ncol(data_grs)) {
    stop("Variable of interest is not the same size as number of observations.")
  }

  # ===========================================================================
  # modify arguments to TMLE functions based on type of type requested
  # ===========================================================================
  if (tmle_type == "glm") {
    if (vim %in% c("ate", "rr")) {
      # set GLM libraries for "tmle" package
      tmle_args$g_lib <- c("SL.mean", "SL.glm")
      tmle_args$Q_lib <- c("SL.mean", "SL.glm")
    }
  }

  # ===========================================================================
  # if NPVI parameter requested, set sensible defaults
  # ===========================================================================
  if (vim == "npvi" & is.null(tmle_args$npvi_descr)) {
    npvi_descr_defaults <- list(
      f = identity, iter = 10, cvControl = 2,
      nMax = 30,
      stoppingCriteria = list(
        mic = 0.001,
        div = 0.001,
        psi = 0.01
      )
    )
    tmle_args$npvi_descr <- npvi_descr_defaults
  }

  # ===========================================================================
  # invoke S4 class constructor for object
  # ===========================================================================
  methy_tmle <- .methytmle(data_grs)
  methy_tmle@call <- call
  methy_tmle@var_int <- var_int

  # ===========================================================================
  # set up parallelization if so desired
  # ===========================================================================
  set_parallel(
    parallel = parallel,
    future_param = future_param,
    bppar_type = bppar_type
  )

  # ===========================================================================
  # screen sites to produce a subset on which to estimate VIMs
  # ===========================================================================
  if (filter == "limma") {
    methy_tmle <- limma_screen(
      methytmle = methy_tmle,
      var_int = var_int,
      type = type
    )
  }

  # ===========================================================================
  # cluster sites based on genomic windows
  # ===========================================================================
  methy_tmle <- cluster_sites(methytmle = methy_tmle)

  # ===========================================================================
  # TMLEs for the Average Treatment Effect (ATE) and Risk Ratio (RR) parameters
  # ===========================================================================
  if (vim %in% c("ate", "rr")) {
    # make sure that the outcome data is of class numeric
    var_of_interest <- as.numeric(methy_tmle@var_int)

    # get names of sites to be added to output object
    cpg_screened_names <- names(methy_tmle[methy_tmle@screen_ind])

    # object of screened CpG site indices to loop over in TMLE procedure
    methy_tmle_ind <- seq_along(methy_tmle@screen_ind)

    # If specifying estimation over a limited number of sites
    if (!is.null(sites_comp)) {
      cpg_screened_names <- cpg_screened_names[seq_len(sites_comp)]
      methy_tmle_ind <- methy_tmle_ind[seq_len(sites_comp)]
    }

    # set the estimation function based on the choice of backend package
    if (tmle_backend == "tmle") {
      methyvim_est <- methyvim_tmle
    } else if (tmle_backend == "drtmle") {
      methyvim_est <- methyvim_tmle
      # methyvim_est <- methyvim_drtmle
    } else if (tmle_backend == "tmle.npvi") {
      methyvim_est <- methyvim_tmle
      # methyvim_est <- methyvim_npvi
    }

    # avoid some try-errors by wrapping estimation function in try statement
    methyvim_est_try <- wrap_in_try(methyvim_est)

    # Perform the estimation procedure in parallel
    methy_vim_out <- BiocParallel::bplapply(methy_tmle_ind,
      FUN = methyvim_est_try,
      methytmle_screened = methy_tmle,
      var_of_interest = var_of_interest,
      type = type,
      corr = corr_max,
      obs_per_covar = obs_per_covar,
      target_param = vim,
      g_lib = tmle_args$g_lib,
      Q_lib = tmle_args$Q_lib,
      cv_folds = tmle_args$cv_folds,
      return_ic = return_ic
    )
    methy_vim_out <- do.call(rbind.data.frame, methy_vim_out)

    # TMLE procedure is now done, so let's just clean up the output
    if (vim == "ate") {
      colnames(methy_vim_out) <- c(
        "lwr_ci", "est_ate", "upr_ci",
        "var_ate", "pval", "n_neighbors", "n_neighbors_control",
        "max_cor_neighbors"
      )
      methy_tmle@param <- "Average Treatment Effect"
    } else {
      colnames(methy_vim_out) <- c(
        "lwr_ci", "est_logrr", "upr_ci",
        "var_logrr", "pval", "n_neighbors", "n_neighbors_control",
        "max_cor_neighbors"
      )
      methy_tmle@param <- "Risk Ratio"
    }
    rownames(methy_vim_out) <- cpg_screened_names
    methy_tmle@vim <- methy_vim_out
  }
  return(methy_tmle)
}

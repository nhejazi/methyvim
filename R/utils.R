#' CpG Neighborhoods from Genomic Distance
#'
#' Clustering of CpG sites to define CpG neighborhoods based on distance (bp).
#'
#' @param methytmle Object of class \code{methytmle} produced from an object of
#'        class \code{GenomicRatioSet}, but containing extra slots.
#' @param window_size Numeric giving the number of base pairs used to define
#'        neighborhoods along the genome. CpG sites within this distance (bp) of
#'        one another are denoted as neighbors. Chromosome boundaries and other
#'        biological constraints are respected (see the documentation available
#'        for \code{bumphunter::clusterMaker} for details).
#'
#' @return An object of class \code{methytmle} with the "clusters" slot filled
#'        in. The "clusters" slot contains a \code{numeric} vector as long as
#'        the number of CpG sites. Each entry in the vector is a neighborhood
#'        assignment used in the estimation procedure.
#'
#' @importFrom bumphunter clusterMaker
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom IRanges ranges
#
cluster_sites <- function(methytmle, window_size = 1000) {
  gr <- SummarizedExperiment::rowRanges(methytmle)
  pos <- BiocGenerics::start(IRanges::ranges(gr))
  clusters <- bumphunter::clusterMaker(chr = GenomeInfoDb::seqnames(gr),
                                       pos = pos,
                                       assumeSorted = FALSE,
                                       maxGap = window_size)
  methytmle@clusters <- as.numeric(clusters)
  return(methytmle)
}

################################################################################

#' FDR-MSA correction
#'
#' Modified FDR Controlling Procedure for Multi-Stage Analyses (MJ van der Laan
#' and C Tuglus, 2009, <doi:10.2202/1544-6115.1397>)
#'
#' @param pvals Numeric vector containing the p-values that result from any
#'        chosen statistical hypothesis testing procedure.
#' @param total_obs Numeric indicating the total number of observations that
#'        would have been available for testing prior to the selection procedure
#'        employed in the multi-stage analysis performed.
#'
#' @return A \code{numeric} vector of corrected p-values, controlling the False
#'         Discovery Rate, using the method of Tuglus and van der Laan.
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
#
fdr_msa <- function(pvals, total_obs) {
  pvals_not_tested <- rep(1, total_obs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fdr_adj <- stats::p.adjust(pvals_all, method = "fdr")
  fdr_out <- fdr_adj[seq_along(pvals)]
  return(fdr_out)
}

################################################################################

#' Easily set up parallelization
#'
#' @param parallel Logical indicating whether parallelization ought to be used.
#'        Parallelization is invoked via a combination of \code{BiocParallel}
#'        and \code{future}. If \code{TRUE} the default method uses multiprocess
#'        evaluation, though other \code{future::plan}s may be specified using
#'        an optional argument. If \code{FALSE}, sequential computation is used
#'        and a warning message is issued.
#' @param future_param Character (if not \code{NULL}) specifying a particular
#'        parallelization approach to be used. For a list of the options, see
#'        the documentation for \code{future::plan}. If the previous argument
#'        (\code{parallel}) is set to \code{FALSE}, this argument is ignored and
#'        sequential computation is invoked via \code{future::sequential}.
#' @param bppar_type Character specifying the type of backend to be used with
#'        the parallelization invoked by \code{BiocParallel}. Consult the manual
#'        page for \code{BiocParallel::BiocParallelParam} for possible types and
#'        descriptions on their appropriate uses. The default for this argument
#'        is \code{NULL}, which silently uses \code{BiocParallel::DoparParam}.
#'
#' @return Nothing. This function is designed to be called for its side-effect
#'         of registering a parallel backend (for \code{BiocParallel}) and/or
#'         \code{future::plan}, making parallel computation a trivial process.
#'
#' @importFrom BiocParallel register bpprogressbar DoparParam
#' @importFrom future plan multiprocess sequential
#' @importFrom doFuture registerDoFuture
#
set_parallel <- function(parallel = c(TRUE, FALSE),
                         future_param = NULL,
                         bppar_type = NULL) {
  # invoke a future-based backend
  doFuture::registerDoFuture()

  if (parallel == TRUE) {
    if (!is.null(future_param)) {
      set_future_param <- parse(text = paste0("future", "::", future_param))
      future::plan(eval(set_future_param))
    } else {
      future::plan(future::multiprocess)
    }
  } else if (parallel == FALSE) {
    warning(paste("Sequential evaluation is strongly discouraged.",
                  "\n Proceed with caution."))
    future::plan(future::sequential)
  }
  if (!is.null(bppar_type)) {
    bp_type <- eval(parse(text = paste0("BiocParallel", "::",
                                        bppar_type, "()")))
  } else {
    bp_type <- BiocParallel::DoparParam()
  }
  # try to use a progress bar is supported in the parallelization plan
  BiocParallel::bpprogressbar(bp_type) <- TRUE
  # register the chosen parallelization plan
  BiocParallel::register(bp_type, default = TRUE)
}

################################################################################

#' Enforce the Assumption of Positivity
#'
#' Discretize continuous variables in the adjustment set (W) of a TMLE procedure
#' in order to avoid practical violations of the assumption of positivity.
#' Discretizes the numeric columns of an input matrix such that the newly
#' created levels of each variable individually contain at least a specified
#' mass when considering each level against levels of the treatment variable.
#'
#' @param A Numeric giving the levels of the (discretized) treatment variable.
#' @param W Data.Frame or Matrix containing the covariates in the adjustment set
#'        to be discretized against the levels of the treatment variable.
#' @param pos_min Numeric indicating the minimum mass (as a proportion) of the
#'        observations to be included in any cell of the table composed of the
#'        levels of the treatment against levels of an adjustment covariate.
#' @param q_init Numeric indicating the initial number of levels to discretize a
#'        given adjustment variable into. This defaults to quantiles.
#'
#' @return A numeric vector with the adjustment variables re-coded into discrete
#'         levels respecting the minimum mass requested in each table comparing
#'         levels of the treatment against levels of an adjustment covariate.
#'
#' @importFrom gtools quantcut
#
force_positivity <- function(A, W, pos_min = 0.1, q_init = 10) {
  stopifnot(length(A) == nrow(W))

  if (class(W) != "data.frame") W <- as.data.frame(W) # cover use of "ncol"
  out_w <- NULL # concatenate W columnwise as we discretize each covar below

  for (obs_w in seq_len(ncol(W))) {
    in_w <- as.numeric(W[, obs_w])
    discr_w <- as.numeric(as.factor(gtools::quantcut(x = in_w, q = q_init)))
    check <- sum((table(A, discr_w) / length(A)) < pos_min)
    next_guess_q <- q_init
    while (check > 0) {
      next_guess_q <- (next_guess_q - 1)
      discr_w <- as.numeric(as.factor(gtools::quantcut(x = in_w,
                                                       q = next_guess_q)))
      check <- sum((table(A, discr_w) / length(A)) < pos_min)
    }
    out_w <- cbind(out_w, discr_w)
  }
  out <- as.data.frame(out_w)
  colnames(out) <- colnames(W)
  rownames(out) <- rownames(W)
  if(length(which(colSums(out) == length(A))) > 0) {
    out <- out[, -which(colSums(out) == length(A)), drop = FALSE]
  }
  return(out)
}


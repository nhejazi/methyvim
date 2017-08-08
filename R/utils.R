#' CpG Neighborhoods from Genomic Distance
#'
#' Clustering of CpG sites to define CpG neighborhoods based on distance (bp).
#'
#' @param methy_tmle Object of class \code{methytmle} produced from an object of
#'        class \code{GenomicRatioSet}, but containing extra slots.
#' @param window_size Numeric giving the number of base pairs used to define
#'        neighborhoods along the genome. CpG sites within this distance (bp) of
#'        one another are denoted as neighbors. Chromosome boundaries and other
#'        biological constraints are respected (see the documentation available
#'        for \code{bumphunter::clusterMaker} for details).
#'
#' @importFrom bumphunter clusterMaker
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
#' @importFrom IRanges ranges
#'
cluster_sites <- function(methy_tmle, window_size = 1000) {
  gr <- SummarizedExperiment::rowRanges(methy_tmle)
  pos <- BiocGenerics::start(IRanges::ranges(gr))
  clusters <- bumphunter::clusterMaker(chr = GenomeInfoDb::seqnames(gr),
                                       pos = pos,
                                       assumeSorted = FALSE,
                                       maxGap = window_size)
  methy_tmle@clusters <- as.numeric(clusters)
  return(methy_tmle)
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
#' @importFrom stats p.adjust
#'
fdr_msa <- function(pvals, total_obs) {
  pvals_not_tested <- rep(1, total_obs - length(pvals))
  pvals_all <- c(pvals, pvals_not_tested)
  fdr_adj <- p.adjust(pvals_all, method = "fdr")
  fdr_out <- fdr_adj[seq_len(pvals)]
  return(fdr_out)
}

################################################################################

#' Set up parallelization
#'
#' @param parallel Numeric or Logical indicating either the number of cores to
#'        use in parallelized computation or whether parallelization ought to be
#'        used at all. If \code{TRUE}, all available cores are found via a call
#'        to \code{parallel::detectCores} and used. If \code{FALSE}, computation
#'        is performed sequentially (with single core) and a warning is issued.
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
set_parallel <- function(parallel) {
  if (class(parallel) == "numeric") doParallel::registerDoParallel(parallel)
  if (class(parallel) == "logical") {
    nCores <- parallel::detectCores()
    if (nCores > 1) {
      doParallel::registerDoParallel(nCores)
    } else {
      warning("option 'parallel' is set to TRUE but only 1 core detected.")
    }
    if (parallel == FALSE) {
      warning("parallelization has been set to FALSE: the estimation procedure
               will likely take on the order of days to run to completion.")
    }
  }
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
#'
#' @export
#'
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

#' Clustering of sites to generate CpG neighborhoods
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

#' Enforce the Assumption of Positivity to Compute the Average Treatment Effect
#'
#' Discretizes the numeric columns of an input matrix such that the newly
#' created levels of each variable individually contain at least a specified
#' mass when considering each level against levels of the treatment variable.
#'
#' @param A...
#' @param W...
#' @param pos_min Numeric indicating...
#' @param q_init Numeric indicating...
#'
#' @return A numeric vector with the data re-coded to based on the quantiles.
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

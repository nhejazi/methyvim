# catch inputs
catch_inputs <- list(data = data_grs, var = var_int, cpg_iss = cpg_is,
                     type = type, vim = vim, neighbors = neighbors,
                     normalize = normalize, filter = filter,
                     min_sites = min_sites, family = family,
                     g_lib = g_lib, Q_lib = Q_lib, parallel = parallel,
                     return_ic = return_ic, shrink_ic = shrink_ic)

# define methytmle class
.methytmle <- methods::setClass(
       Class = "methytmle",
       slots = list(call = "call",
                    screen_ind = "numeric",
                    clusters = "numeric",
                    g = "matrix",
                    Q = "matrix",
                    ic = "data.frame",
                    vim = "data.frame"),
       contains = "GenomicRatioSet"
)

# set up methyvim object
call <- "testing123"
class(call) <- "call" # hacking the call object
methy_tmle <- .methytmle(catch_inputs$data)
methy_tmle@call <- call

# using LIMMA for screening
limma_screen <- function(methytmle, var_int, type, cutoff = 0.05) {
  # setup design matrix
  design <- as.numeric(colData(methytmle)[, var_int])
  design <- as.matrix(cbind(rep(1, times = length(design)), design))

  # create expression object for modeling
  if (type == "Beta") {
    methytmle_exprs <- minfi::getBeta(methytmle)
  } else if (type == "Mval") {
    methytmle_exprs <- minfi::getM(methytmle)
  }

  # fit limma model and apply empirical Bayes shrinkage
  mod_fit <- limma::lmFit(object = methytmle_exprs, design = design)
  mod_fit <- limma::eBayes(mod_fit)

  # extract indices of relevant CpG sites
  tt_out <- limma::topTable(mod_fit, coef = 2, num = Inf, sort.by = "none")
  indices_pass <- which(tt_out$P.Value < cutoff)

  # add to appropriate slot in the methytmle input object
  methytmle@screen_ind <- indices_pass
  return(methytmle)
}

# function to cluster sites
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

# force positivity assumption to hold
force_positivity <- function(A, W, pos_min = 0.1) {
  stopifnot(length(A) == nrow(W))

  n_obs <- length(A)
  guess_w <- round((pos_min * n_obs) / length(unique(A))) # heuristic W binning

  if (class(W) != "data.frame") W <- as.data.frame(W) # cover use of "ncol"
  out_w <- NULL # concatenate W columnwise as we discretize each covar below

  for (obs_w in seq_len(ncol(W))) {
    in_w <- as.numeric(W[, obs_w])
    discr_w <- as.numeric(as.factor(gtools::quantcut(x = in_w, q = guess_w)))
    check <- sum((table(A, discr_w) / n_obs) < pos_min)
    next_guess_w <- guess_w
    while (check > 0) {
      next_guess_w <- (next_guess_w - 1)
      discr_w <- as.numeric(as.factor(gtools::quantcut(x = in_w,
                                                       q = next_guess_w)))
      check <- sum((table(A, discr_w) / n_obs) < pos_min)
    }
    out_w <- cbind(out_w, discr_w)
  }
  out <- as.data.frame(out_w)
  colnames(out) <- colnames(W)
  rownames(out) <- rownames(W)
  if(length(which(colSums(out) == n_obs)) > 0) {
    out <- out[, -which(colSums(out) == n_obs), drop = FALSE]
  }
  return(out)
}

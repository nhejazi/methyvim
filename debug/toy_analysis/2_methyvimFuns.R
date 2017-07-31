# methyvim args
var_int = 3 #exposed
cpg_is = "exposure"
neighbors = 1e3
normalize = NULL
filter = TRUE
min_sites = 1e4
family = "gaussian"
g_lib = c("SL.mean", "SL.glm", "SL.randomForest")
Q_lib = c("SL.mean", "SL.randomForest")
parallel = TRUE
return_ic = TRUE
shrink_ic = TRUE
type = "M"
vim = "npvi"
data_grs <- buccal_funnorm

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

# screen CpG sites using LIMMA method
methy_tmle_screened <- limma_screen(methytmle = methy_tmle,
                                    var_int = catch_inputs$var,
                                    type = catch_inputs$type)

# NOTE: work out ATE procedure "by hand"

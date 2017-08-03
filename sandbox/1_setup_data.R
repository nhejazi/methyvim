# needed pkgs
library(here)
library(minfi)
library(parallel)
library(doParallel)
n_cores <- detectCores()
registerDoParallel(n_cores)

# check system
stopifnot(
  stringr::str_split(Sys.info()["nodename"], "\\.")[[1]][1] == "bluevelvet"
)

# load sample data for playing around
load(here("..", "data", "tmleMethyData", "epigen_grs.RData"))

# methyvim args
exp_var_int = 1 # binary exposure (for ATE)
out_var_int = 2 # continuous outcome (for NPVI)
window_bp = 1e3
corr_max = 0.35
obs_per_covar = 15
type = "Mval"
filter = "limma"
cutoff = 0.05
parallel = TRUE
return_ic = FALSE
shrink_ic = FALSE
family = "binomial"  # for logistic fluctuation model with scaled Y
g_lib = c("SL.mean", "SL.glm", "SL.randomForest")
Q_lib = c("SL.mean", "SL.randomForest")
npvi_cutoff = 0.25
npvi_descr = list(f = identity, iter = 10, cvControl = 2, nMax = 30,
                  stoppingCriteria = list(mic = 0.001, div = 0.001,
                                          psi = 0.01)
                 )

# catch inputs for ATE procedure
catch_inputs_ate <- list(data = grs, var = exp_var_int, vim = "ATE",
                         type = type, window = window_bp, corr = corr_max,
                         obs_per_covar = obs_per_covar, par = parallel,
                         filter = filter, filter_cutoff = cutoff,
                         return_ic = return_ic, shrink_ic = shrink_ic,
                         tmle_args = list(family = family,
                                          g_lib = g_lib,
                                          Q_lib = Q_lib,
                                          npvi_cutoff = npvi_cutoff,
                                          npvi_descr = npvi_descr
                                         )
                        )

# catch inputs for NPVI procedure
catch_inputs_npvi <- list(data = grs, var = exp_var_int, vim = "NPVI",
                          type = type, window = window_bp, corr = corr_max,
                          obs_per_covar = obs_per_covar, par = parallel,
                          filter = filter, filter_cutoff = cutoff,
                          return_ic = return_ic, shrink_ic = shrink_ic,
                          tmle_args = list(family = family,
                                           g_lib = g_lib,
                                           Q_lib = Q_lib,
                                           npvi_cutoff = npvi_cutoff,
                                           npvi_descr = npvi_descr
                                          )
                         )

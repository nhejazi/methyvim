# needed pkgs
library(here)
library(minfi)
library(foreach)
library(parallel)
library(doParallel)
n_cores <- detectCores()
registerDoParallel(n_cores)

# check system, makes sure you're on bluevelvet
stopifnot(
  stringr::str_split(Sys.info()["nodename"], "\\.")[[1]][1] == "bluevelvet"
)

# load sample data for playing around
load(here("..", "data", "tmleMethyData", "epigen_grs.RData"))

# methyvim args, defaults that are inputted into the main functions
exp_var_int = 1 # binary exposure (for ATE)
# column variables in the coldata matrix
out_var_int = 2 # continuous outcome (for NPVI)
window_bp = 1e3
corr_max = 0.35
# maximum correlation that a neighboring site can have with the target site
obs_per_covar = 15
# how many obs you need for every covariate you include in W (20/covariate)
# otherwise you have insufficient data to control for all of this
type = "Mval"
filter = "limma"
# one based on npvi and another for wilsons package, pass object into screening
# appropriate function and sites that pass are the ones we analyze
cutoff = 0.05 # defines what passes through the filter, p-value
# limma might want to use t-stat
# npvi use p-value
# data adaptive w wilon might want to use rank of site w/i final table produced
parallel = TRUE
return_ic = FALSE # ic estimates out, slot in methytmle object
shrink_ic = FALSE # apply limma to reduce variance in ic based estimates
# at that point you wont care about values from tmle estimation procedure
# ic based estimate for each site and variance of these sites is variance of
# efficient influence curve, so each site is distributed w mean (psi)
# then do t-test w the expectation with the ate
# (controlling erratic spikes in the t stat based on low variance)
family = "binomial"  # for logistic fluctuation model with scaled Y
# or gaussian
g_lib = c("SL.mean", "SL.glm", "SL.randomForest") #sensible defaults
Q_lib = c("SL.mean", "SL.randomForest")
npvi_cutoff = 0.25
# keep a continuous (ate discretizes things) but you have to map some range of
# the values to zero (null range with enough observations)
# cutoff defines the quantile in the obs data that you'll use to map to zero
# antione had long discussion about this and theres an optimal way to deciding
# this cutoff
npvi_conf = 0.95
npvi_descr = list(f = identity, iter = 10, cvControl = 2, nMax = 30,
                  stoppingCriteria = list(mic = 0.001, div = 0.001,
                                          psi = 0.01)
                 )
# below is done automattically, just done for organization
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
catch_inputs_npvi <- list(data = grs, var = out_var_int, vim = "NPVI",
                          type = type, window = window_bp, corr = corr_max,
                          obs_per_covar = obs_per_covar, par = parallel,
                          filter = filter, filter_cutoff = cutoff,
                          return_ic = return_ic, shrink_ic = shrink_ic,
                          tmle_args = list(family = family,
                                           g_lib = g_lib,
                                           Q_lib = Q_lib,
                                           npvi_cutoff = npvi_cutoff,
                                           npvi_conf = npvi_conf,
                                           npvi_descr = npvi_descr
                                          )
                         )

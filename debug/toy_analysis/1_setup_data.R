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
load(here("..", "data", "epic-smith-2017",
          "buccal_epic_funnorm_methyvim.RData"))

# methyvim args
var_int = 3 #exposed
window_bp = 1e3
corr_max = 0.35
cutoff = 0.05
obs_per_var = 15
preprocess = NULL
filter = TRUE
family = "binomial"  # for logistic fluctuation model with scaled Y
g_lib = c("SL.mean", "SL.glm", "SL.randomForest")
Q_lib = c("SL.mean", "SL.randomForest")
parallel = TRUE
return_ic = TRUE
shrink_ic = TRUE
dimen_red = TRUE
type = "Mval"
vim = "npvi"
data_grs <- buccal_funnorm

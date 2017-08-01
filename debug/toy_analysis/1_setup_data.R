# needed pkgs
library(here)
library(minfi)

# check system
stopifnot(
  stringr::str_split(Sys.info()["nodename"], "\\.")[[1]][1] == "bluevelvet"
)

# load sample data for playing around
load(here("..", "data", "epic-smith-2017",
          "buccal_epic_funnorm_methyvim.RData"))

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

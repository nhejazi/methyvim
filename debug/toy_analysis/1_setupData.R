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

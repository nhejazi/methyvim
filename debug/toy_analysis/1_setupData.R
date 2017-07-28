# needed pkgs
library(here)
library(minfi)

# check system
stop(Sys.info()["nodename"] != "bluevelvet.biostat.berkeley.edu")

# load sample data for playing around
load(here("..", "data", "epic-smith-2017",
          "buccal_epic_funnorm_methyvim.RData"))

# ...

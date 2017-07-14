# setup project
library(here)

# setting project and data directories...
proj_dir <- here()
if (substr(Sys.info()["nodename"], 1, 10) == "bluevelvet") {
  data_mtmle <- normalizePath(here("..", "..", "..", "..",
                                   "local_ATE-methylation",
                                   "methytmle_data_fixed.RData"))
} else {
  warning("Use Bluevelvet server to analyze this data set (due to size).")
}
load(data_mtmle)


# extract annotation and M-values
annot <- methytmleTestData$annot
batch <- methytmleTestData$batch
mvals <- methytmleTestData$mvals[, as.character(batch$Sample_Name)]

# extract outcomes of interest
outcomes <- c("fsiq_i_7y", "vciq_i_7y", "priq_i_7y", "wmiq_i_7y", "psiq_7y",
              "baattss_7y", "pbdss_7y", "ftdr_7y")

# extract potential confounders
confounders <- c("educcat", "sppstd_6m", "CD8T.bakulski", "CD4T.bakulski",
                 "NK.bakulski", "Bcell.bakulski", "Mono.bakulski",
                 "Gran.bakulski", "nRBC.bakulski")

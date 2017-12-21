.onAttach <- function(...) {
  packageStartupMessage(paste("methyvim: Targeted Data-Adaptive Estimation and",
                              "Inference",
                              "\n for Differential Methylation Analysis"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}

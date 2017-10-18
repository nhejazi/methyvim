.onAttach <- function(...) {
  packageStartupMessage(paste("methyvim: Targeted Learning for Differential",
                              "Methylation Analysis"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}

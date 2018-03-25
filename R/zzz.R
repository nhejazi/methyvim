.onAttach <- function(...) {
  packageStartupMessage(paste0("methyvim v",
                               utils::packageDescription("methyvim")$Version,
                               ": Targeted Variable Importance for ",
                               "Differential Methylation Analysis"))
}

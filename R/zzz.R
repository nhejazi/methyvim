.onAttach <- function(...) {
  packageStartupMessage("methyvim: Nonparametric Variable Importance for Differential Methylation Analysis")
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}

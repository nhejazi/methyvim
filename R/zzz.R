.onAttach <- function(...) {
  packageStartupMessage("methyvim: Differential Methylation Analysis with Nonparametric Variable Importance Measures")
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}

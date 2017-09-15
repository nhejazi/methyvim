.onAttach <- function(...) {
  packageStartupMessage("methyvim: Nonparametric Differential Methylation Analysis with Variable Importance Measures")
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}

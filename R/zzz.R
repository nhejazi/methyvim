.onAttach <- function(...) {
  packageStartupMessage(paste("methyvim: Nonparametric Differential",
                              "Methylation Analysis \n with Targeted Estimates",
                              "of Variable Importance Measures"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("methyvim")$Version)
}


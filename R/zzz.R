.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "methyvim v",
    utils::packageDescription("methyvim")$Version,
    ": Targeted, Robust, and Model-free Differential Methylation Analysis"
  ))
}

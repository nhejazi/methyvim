.onAttach <- function(...) {
    packageStartupMessage("methadapt: targeted learning for DNA methylation")
  packageStartupMessage("Version: ",      utils::packageDescription("methadapt")$Version)
}

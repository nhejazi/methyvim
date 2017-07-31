do_preprocess <- function(methy_tmle, type = "Funnorm") {
  #norm_type <- paste0("minfi", "::", "preprocess", type, ...)
  #methy_tmle_norm <- eval(parse(text = paste(norm_type, ""))
  methy_tmle_norm <- minfi::preprocessFunnorm(methy_tmle, bgCorr = TRUE,
                                              dyeCorr = TRUE,
                                              ratioConvert = FALSE)
  methy_tmle_norm <- minfi::ratioConvert(methy_tmle, what = "both",
                                         keepCN = TRUE)
  # map to genome
  methy_tmle_norm_mapped <- minfi::mapToGenome(methy_tmle_norm)

  # drop SNPs
  methy_tmle_out <- minfi::addSnpInfo(methy_tmle_norm_mapped)
  methy_tmle_out <- minfi::dropLociWithSnps(methy_tmle_out,
                                            snps = c("SBE", "CpG"),
                                            maf = 0)
  return(methy_tmle_out)
}

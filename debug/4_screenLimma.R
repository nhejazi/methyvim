# use LIMMA to assess neurocognitive scores using fully adjusted models
tt_neuro_full <- vector("list", length(outcomes))

for (neuro in 1:length(outcomes)) {
  # normalization and subject reduction
  meth <- as.data.frame(mcols(cpg_gr_mval)[, 2:length(mcols(cpg_gr_mval))])
  meth_complete <- meth[, keep_cases[[neuro]]]

  # fit LIMMA model
  fit_full <- limma::lmFit(meth_complete, designMats_full[[neuro]])
  fit_full <- limma::eBayes(fit_full)
  tt_full <- limma::topTable(fit_full, coef = 2, num = Inf, sort.by = "none")

  # store topTable result for later examination
  tt_neuro_full[[neuro]] <- tt_full
}

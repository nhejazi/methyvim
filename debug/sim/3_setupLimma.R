# use LIMMA to assess association with neurodevelopment with adjusted models
designMats_full <- vector("list", length(outcomes))
designMats_min <- vector("list", length(outcomes))
keep_cases <- vector("list", length(outcomes))

for (covar in 1:length(outcomes)) {
  # extract full design matrix
  covarsNames <- c(outcomes[covar], confounders, "batch")
  design <- subset(batch, select = covarsNames)

  # remove cases where covariate of interest is missing
  keep_cases[[covar]] <- complete.cases(design)
  design <- design[keep_cases[[covar]], ]

 # set up intercept term
 design <- as.data.frame(cbind(rep(1, nrow(design)), design))
 colnames(design)[1] <- "intercept"
 design <- sapply(design, as.numeric)

 # reduce to minimal design matrix
 design_min <- subset(design, select = c("intercept", outcomes[covar]))

 # store design matrix for later use
 designMats_full[[covar]] <- design
 designMats_min[[covar]] <- design_min
}

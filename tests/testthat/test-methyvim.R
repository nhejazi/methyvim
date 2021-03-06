context("Wrapper function for TMLE for supported target parameters")

# libraries and data
suppressMessages(library(tmle))
suppressMessages(library(minfi))
suppressMessages(library(SummarizedExperiment))
library(methyvimData)
data(grsExample)
var_int <- colData(grsExample)[, 1]


# run TMLE procedure for the ATE parameter over M-values with Limma filtering
methyvim_out_ate <- suppressWarnings(
  methyvim(
    data_grs = grsExample, var_int = var_int, vim = "ate",
    type = "Mval", filter = "limma", filter_cutoff = 0.05,
    sites_comp = 3, parallel = FALSE, tmle_type = "sl"
  )
)

# run TMLE procedure for the RR parameter over M-values with Limma filtering
methyvim_out_rr <- suppressWarnings(
  methyvim(
    data_grs = grsExample, var_int = var_int, vim = "rr",
    type = "Mval", filter = "limma", filter_cutoff = 0.05,
    sites_comp = 3, parallel = FALSE, tmle_type = "sl"
  )
)


# test that output object is of the appropriate type and class
test_that("Output object is of class type S4", {
  expect_equal(typeof(methyvim_out_ate), "S4")
})

test_that("Output object is of class methytmle", {
  expect_equivalent(class(methyvim_out_ate), "methytmle")
})

test_that("Output object inherits from minfi::GenomicRatioSet", {
  expect_is(methyvim_out_ate, "GenomicRatioSet")
})


# test that output object contains results of the appropriate form
test_that("Slot call of methytmle object contains object of class call", {
  expect_equal(class(methyvim_out_ate@call), "call")
})

test_that("Slot parameter lists the ATE as the target parameter", {
  expect_equal(methyvim_out_ate@param, "Average Treatment Effect")
})

test_that("Slot of screening IDs are of numeric type and correct length", {
  expect_equal(class(methyvim_out_ate@screen_ind), "numeric")
  expect_equal(length(methyvim_out_ate@screen_ind), 28)
})

test_that("Variable importance results in tolerance and of correct length", {
  expect_equal(nrow(methyvim_out_ate@vim), 3)
  expect_lt(
    sum(range(methyvim_out_ate@vim$pval) - c(0.225938, 0.470252)),
    0.02
  )
})

test_that("Cluster IDs are of appropriate length and unique length", {
  expect_equal(length(methyvim_out_ate@clusters), nrow(methyvim_out_ate))
  expect_lt(length(unique(methyvim_out_ate@clusters)), nrow(methyvim_out_ate))
})

test_that("Variable importance results within tolerance for the RR parameter", {
  expect_lt(
    sum(range(methyvim_out_rr@vim$pval) - c(0.2281579, 0.4707649)),
    0.02
  )
})

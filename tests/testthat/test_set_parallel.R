context("parallelization with futures and BiocParallel")
suppressMessages(library(future))
suppressMessages(library(doFuture))
suppressMessages(library(BiocParallel))


# test: BiocParallel::DoparParam is invoked under parallel=TRUE
set_parallel(parallel = TRUE)
test_that("registers BiocParallel::DoparParam by default for parallel=TRUE", {
  expect_equivalent(names(registered()), "DoparParam")
})

test_that("registers future::multiprocess by default for parallel=TRUE", {
  expect_true("multiprocess" %in% class(plan()))
})


# test: registers a future plan "sequential" with parallel=FALSE
suppressWarnings(set_parallel(parallel = FALSE))
test_that("registers future::sequential for parallel=FALSE", {
  expect_true("sequential" %in% class(plan()))
})


# test: registers backends similarly to BiocParallel
set_parallel(parallel = TRUE, bppar_type = "SnowParam")
test_that("registers BiocParallel::SnowParam when so asked", {
  expect_true("SnowParam" %in% names(registered()))
})

set_parallel(parallel = TRUE, bppar_type = "MulticoreParam")
test_that("registers BiocParallel::MulticoreParam when so asked", {
  expect_true("MulticoreParam" %in% names(registered()))
})

set_parallel(parallel = TRUE, bppar_type = "BatchJobsParam")
test_that("registers BiocParallel::BatchJobsParam when so asked", {
  expect_true("BatchJobsParam" %in% names(registered()))
})

set_parallel(parallel = TRUE, bppar_type = "SerialParam")
test_that("registers BiocParallel::SerialParam when so asked", {
  expect_true("SerialParam" %in% names(registered()))
})


# test: registers future plans exactly as future::plan does
set_parallel(parallel = TRUE, future_param = "sequential")
test_that("sets up sequential future plan when so asked", {
  expect_true("sequential" %in% class(plan()))
})

set_parallel(parallel = TRUE, future_param = "transparent")
test_that("sets up transparent future plan when so asked", {
  expect_true("transparent" %in% class(plan()))
})


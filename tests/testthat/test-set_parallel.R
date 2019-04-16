context("parallelization with futures and BiocParallel")
library(future)
library(doFuture)
library(BiocParallel)


# test: BiocParallel::DoparParam is invoked under parallel=TRUE
test_that("registers BiocParallel::DoparParam by default for parallel=TRUE", {
  methyvim:::set_parallel(parallel = TRUE)
  expect_true("DoparParam" %in% names(registered()))
})

test_that("registers future::multiprocess by default for parallel=TRUE", {
  methyvim:::set_parallel(parallel = TRUE)
  expect_true("multiprocess" %in% class(plan()))
})


# test: registers a future plan "sequential" with parallel=FALSE
test_that("registers future::sequential for parallel=FALSE", {
  suppressWarnings(methyvim:::set_parallel(parallel = FALSE))
  expect_true("sequential" %in% class(plan()))
})


# test: registers backends similarly to BiocParallel
test_that("registers BiocParallel::SnowParam when so asked", {
  methyvim:::set_parallel(parallel = TRUE, bppar_type = "SnowParam")
  expect_true("SnowParam" %in% names(registered()))
})

test_that("registers BiocParallel::MulticoreParam when so asked", {
  skip_on_os("windows")
  methyvim:::set_parallel(parallel = TRUE, bppar_type = "MulticoreParam")
  expect_true("MulticoreParam" %in% names(registered()))
})

test_that("registers BiocParallel::BatchJobsParam when so asked", {
  methyvim:::set_parallel(parallel = TRUE, bppar_type = "BatchJobsParam")
  expect_true("BatchJobsParam" %in% names(registered()))
})

test_that("registers BiocParallel::SerialParam when so asked", {
  skip_on_os("windows")
  methyvim:::set_parallel(parallel = TRUE, bppar_type = "SerialParam")
  expect_true("SerialParam" %in% names(registered()))
})


# test: registers future plans exactly as future::plan does
test_that("sets up sequential future plan when so asked", {
  methyvim:::set_parallel(parallel = TRUE, future_param = "sequential")
  expect_true("sequential" %in% class(plan()))
})

test_that("sets up transparent future plan when so asked", {
  methyvim:::set_parallel(parallel = TRUE, future_param = "transparent")
  expect_true("transparent" %in% class(plan()))
})

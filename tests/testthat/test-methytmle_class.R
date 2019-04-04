context("methyvim class contains expected properties (subclasses, methods)")

# construct methytmle object from methyvimData example data
suppressMessages(library(minfi))
library(methyvimData)
data(grsExample)
methy_tmle <- .methytmle(grsExample)


# test: is of class methytmle and type S4
test_that("methy_tmle object created is of class methytmle", {
  expect_equivalent(class(methy_tmle), "methytmle")
})

test_that("methytmle objects are of the S4 general class type", {
  expect_equivalent(typeof(methy_tmle), "S4")
})

test_that("methytmle objects inherit from minfi::GenomicRatioSet", {
  expect_is(methy_tmle, "GenomicRatioSet")
})


# test: subclass inherits from minfi objects
test_that("methytmle objects inherit from minfi::GenomicRatioSet", {
  expect_true(inherits(methy_tmle, "GenomicRatioSet"))
})


# test: method "param" extracts an appropriate character (empty) for an
#       initialized methytmle object (without running estimation)
test_that("accessor function 'param' for methytmle class acts appropriately", {
  param_name <- "Just Another Parameter"
  methy_tmle@param <- param_name
  expect_identical(param(methy_tmle), param_name)
})

# test: method "vim" extracts an appropriate data.frame (i.e., empty) for an
#       initialized methytmle object (without running estimation)
test_that("accessor function 'vim' for methytmle class acts appropriately", {
  example_vim <- as.data.frame(replicate(3, rnorm(nrow(grsExample))))
  methy_tmle@vim <- example_vim
  expect_identical(vim(methy_tmle), example_vim)
})

context("methyvim class contains expected properties (subclasses, etc.)")

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

context("screening procedure based on limma")

# use example data from methyvimData package
library(methyvimData)
data(grsExample)
var_int <- colData(grsExample)[, 1]
methytmle <- .methytmle(grsExample)
screened <- methyvim:::limma_screen(methytmle, var_int = var_int,
                                    type = "Mval")
screened_betas <- methyvim:::limma_screen(methytmle, var_int = var_int,
                                          type = "Beta")


# check the screen_ind slot
test_that("no sites passing screening before running Limma procedure", {
  expect_equal(length(methytmle@screen_ind), 0)
})

test_that("some sites passing screening after running Limma on M-values", {
  expect_gt(length(screened@screen_ind), 0)
})

test_that("some sites passing screening after running Limma on Beta-values", {
  expect_gt(length(screened_betas@screen_ind), 0)
})


# screening does not change the class properties
test_that("screened object is of S4 class system", {
  expect_equivalent(typeof(screened), "S4")
})

test_that("screened object is of class methyvim", {
  expect_equivalent(class(screened), "methytmle")
})

test_that("classes of screened object and of unscreened object match", {
  expect_equivalent(class(screened), class(methytmle))
})

test_that("class types of screened object and of unscreened object match", {
  expect_equivalent(typeof(screened), typeof(methytmle))
})


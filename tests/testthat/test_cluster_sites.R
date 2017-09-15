context("clustering sites using bumphunter")

# use example data from methyvimData package
library(methyvimData)
data(grsExample)
methytmle <- .methytmle(grsExample)
clustered <- methyvim:::cluster_sites(methytmle)


# check cluster IDs
test_that("each site is assigned a cluster ID", {
  expect_equal(length(clustered@clusters), nrow(clustered))
})

test_that("non-unique cluster IDs: fewer clusters than individual sites", {
  expect_lt(length(unique(clustered@clusters)), nrow(clustered))
})


# clustering does not change the class properties
test_that("clustered object is of S4 class system", {
  expect_equivalent(typeof(clustered), "S4")
})

test_that("clustered object is of class methyvim", {
  expect_equivalent(class(clustered), "methytmle")
})

test_that("classes of clustered object and of unclustered object match", {
  expect_equivalent(class(clustered), class(methytmle))
})

test_that("class types of clustered object and of unclustered object match", {
  expect_equivalent(typeof(clustered), typeof(methytmle))
})


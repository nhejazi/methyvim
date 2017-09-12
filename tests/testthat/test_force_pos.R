context("discretization enforces practical positivity assumptions")

# simulate toy data
n <- 1e2
vars <- 5
w <- replicate(vars, rnorm(n))
a <- as.numeric(rowSums(w) > mean(rowSums(w)))

# enforce positivity
discrete_w <- force_positivity(A = a, W = w)

# test: the discretization is sufficiently coarse
test_that("discretized W variable has fewer levels than desired maximum", {
  expect_lt(length(unique(discrete_w[, 1])), 10)
})

# test: check that minimum mass is greater than desired minimum
test_that("minimum mass in discretized W is higher than desired minimum", {
  expect_gte(min(table(discrete_w[, 1], a) / n), 0.10)
})


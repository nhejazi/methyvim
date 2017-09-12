context("FDR-MSA correction")

# simulate toy data
g <- 1e4
n <- 1e2
tests <- 1e2
samp <- sample(seq_len(g), g / tests)
w <- replicate(n, rnorm(g))
a <- cbind(rep(1, n), as.numeric(colSums(w) > mean(colSums(w))))

# fit linear models to the simulated data over a subsample of genes
mod <- lmFit(object = w[samp, ], design = a)
mod <- eBayes(mod)
pvals <- topTable(mod, coef = 2, number = Inf)$P.Value

# test: FDR-MSA with p-values in a "normal" range pushes all p-values to 1
fdrmsa <- fdr_msa(pvals = pvals, total_obs = g)
test_that("corrected p-values are uniquely 1 when not sufficiently low", {
  expect_equal(unique(fdrmsa), 1)
})

# test: number of p-values returned should be same as input
test_that("lengths of input and output match for FDR-MSA", {
  expect_equal(length(fdrmsa), length(pvals))
})

# test: do p-values that are really low get shifted reasonably?
p <- abs(rnorm(n, mean = 1e-8, sd = 1e-2))
fdr_p <- fdr_msa(pvals = p, total_obs = g)

test_that("maximum of corrected p-values is 1", {
  expect_equal(max(range(fdr_p)), 1)
})

test_that("minimum of corrected p-values is less than 1", {
  expect_lt(min(range(fdr_p)), 1)
})


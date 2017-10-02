context("targeting parameters with the classic TMLE functions")

# use example data from methyvimData package
suppressMessages(library(SummarizedExperiment))
library(methyvimData)
data(grsExample)
methytmle <- .methytmle(grsExample)
var_int <- colData(methytmle)[, 1]
screened <- methyvim:::limma_screen(methytmle, var_int = var_int, type = "Mval")
clustered <- methyvim:::cluster_sites(screened)


# consistentcy of output for ATE and RR over target site with M-values
methyvim_mvals_ate <- methyvim:::methyvim_tmle(target_site = 10,
                                               methytmle_screened = clustered,
                                               var_of_interest = var_int,
                                               type = "Mval",
                                               corr = 0.75,
                                               obs_per_covar = 20,
                                               target_param = "ate",
                                               family = "binomial",
                                               return_ic = FALSE
                                              )

test_that("ATE procedure with M-values is consistent for target site", {
  expect_equal(methyvim_mvals_ate,
               c(-0.38903184, -0.17995147, 0.02912890, 0.01137927, 0.09161596,
                 0.000000000, 0.000000000, NA))
})

methyvim_mvals_rr <- methyvim:::methyvim_tmle(target_site = 10,
                                              methytmle_screened = clustered,
                                              var_of_interest = var_int,
                                              type = "Mval",
                                              corr = 0.75,
                                              obs_per_covar = 20,
                                              target_param = "rr",
                                              family = "binomial",
                                              return_ic = FALSE
                                             )

test_that("RR procedure with M-values is consistent for target site", {
  expect_equal(methyvim_mvals_rr,
               c(-0.159795138, -0.073782061, 0.012231016, 0.001925825,
                 0.092706794, 0.000000000, 0.000000000, NA))
})

# consistentcy of output for ATE and RR over target site with Beta-values
methyvim_betas_ate <- methyvim:::methyvim_tmle(target_site = 10,
                                               methytmle_screened = clustered,
                                               var_of_interest = var_int,
                                               type = "Beta",
                                               corr = 0.75,
                                               obs_per_covar = 20,
                                               target_param = "ate",
                                               family = "binomial",
                                               return_ic = FALSE
                                              )

test_that("ATE procedure with Beta-values is consistent for target site", {
  expect_equal(methyvim_betas_ate,
               c(-0.0541869168, -0.0264814441, 0.0012240287, 0.0001998108,
                 0.0610121857, 0.0000000000, 0.0000000000, NA))
})

methyvim_betas_rr <- methyvim:::methyvim_tmle(target_site = 10,
                                              methytmle_screened = clustered,
                                              var_of_interest = var_int,
                                              type = "Beta",
                                              corr = 0.75,
                                              obs_per_covar = 20,
                                              target_param = "rr",
                                              family = "binomial",
                                              return_ic = FALSE
                                             )

test_that("RR procedure with Beta-values is consistent for target site", {
  expect_equal(methyvim_betas_rr,
               c(-0.146258893, -0.071235003, 0.003788888, 0.001465167,
                 0.062742021, 0.000000000, 0.000000000, NA))
})


# simple test of estimation results but at a point with neighbors
#methyvim_neighbors <- methyvim:::methyvim_tmle(target_site = 18,
                                               #methytmle_screened = clustered,
                                               #var_of_interest = var_int,
                                               #type = "Mval",
                                               #corr = 0.75,
                                               #obs_per_covar = 20,
                                               #target_param = "ate",
                                               #family = "binomial",
                                               #return_ic = FALSE
                                              #)

#test_that("Procedure is generally consistent for a site with neighbors", {
  #expect_lt(sum(methyvim_neighbors -
                #c(-0.53382020341, -0.35780299063, -0.18178577786, 0.00806488421,
                  #0.00006769798, 8.0000000000, 8.00000000000, -0.09301901866)),
            #1e-11)
#})


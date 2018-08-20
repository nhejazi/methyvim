## ----reqs, echo=FALSE, message=FALSE---------------------------------------
suppressMessages(library(tmle))
suppressMessages(library(minfi))
suppressMessages(library(SummarizedExperiment))

## ----prelims---------------------------------------------------------------
library(methyvim)
library(methyvimData)

## ----get-data--------------------------------------------------------------
set.seed(479253)
data(grsExample)
grsExample
var_int <- as.numeric(colData(grsExample)[, 1])
table(var_int)

## ----make-methytmle--------------------------------------------------------
mtmle <- .methytmle(grsExample)

## ----methyvim-ate-sl-------------------------------------------------------
suppressMessages(
  methyvim_ate_sl <- methyvim(data_grs = grsExample, sites_comp = 25,
                              var_int = var_int, vim = "ate", type = "Mval",
                              filter = "limma", filter_cutoff = 0.10,
                              parallel = FALSE, tmle_type = "sl"
                             )
)

## ----methyvim-ate-sl-print-------------------------------------------------
vim(methyvim_ate_sl)

## ----methyvim-ate-glm------------------------------------------------------
suppressMessages(
  methyvim_ate_glm <- methyvim(data_grs = grsExample, sites_comp = 25,
                               var_int = var_int, vim = "ate", type = "Mval",
                               filter = "limma", filter_cutoff = 0.10,
                               parallel = FALSE, tmle_type = "glm"
                              )
)

## ----methyvim-ate-glm-print------------------------------------------------
vim(methyvim_ate_glm)

## ----methyvim-rr-sl--------------------------------------------------------
methyvim_rr_sl <- methyvim(data_grs = grsExample, sites_comp = 25,
                            var_int = var_int, vim = "rr", type = "Mval",
                            filter = "limma", filter_cutoff = 0.10,
                            parallel = FALSE, tmle_type = "sl"
                           )

## ----methyvim-rr-sl-print--------------------------------------------------
vim(methyvim_rr_sl)

## ----methyvim-rr-glm-------------------------------------------------------
methyvim_rr_glm <- methyvim(data_grs = grsExample, sites_comp = 25,
                            var_int = var_int, vim = "rr", type = "Mval",
                            filter = "limma", filter_cutoff = 0.10,
                            parallel = FALSE, tmle_type = "glm"
                           )

## ----methyvim-rr-glm-print-------------------------------------------------
vim(methyvim_rr_glm)

## ----setup-minfidata-------------------------------------------------------
suppressMessages(library(minfiData))
data(MsetEx)
mset <- mapToGenome(MsetEx)
grs <- ratioConvert(mset)
grs

## ---- minfidata-maketx-----------------------------------------------------
var_int <- (as.numeric(as.factor(colData(grs)$status)) - 1)
table(var_int)

## ---- minfidata-methyvim---------------------------------------------------
suppressMessages(
  methyvim_cancer_ate <- methyvim(data_grs = grs, var_int = var_int,
                                  vim = "ate", type = "Beta", filter = "limma",
                                  filter_cutoff = 0.20, obs_per_covar = 2,
                                  parallel = FALSE, sites_comp = 125,
                                  tmle_type = "glm"
                                 )
)

## ----vim-cancer-ate--------------------------------------------------------
methyvim_cancer_ate

## ----fdr-msa---------------------------------------------------------------
fdr_p <- fdr_msa(pvals = vim(methyvim_cancer_ate)$pval,
                 total_obs = nrow(methyvim_cancer_ate))

## ----methyvim-pvals-both---------------------------------------------------
plot(methyvim_cancer_ate)

## ----methyvim-volcano------------------------------------------------------
methyvolc(methyvim_cancer_ate)

## ----methyvim-heatmap------------------------------------------------------
methyheat(methyvim_cancer_ate)

## ----session-info, echo=FALSE----------------------------------------------
sessionInfo()


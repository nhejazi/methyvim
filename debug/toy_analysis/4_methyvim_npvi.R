# screen CpG sites using LIMMA method
methy_tmle_screened <- limma_screen(methytmle = methy_tmle,
                                    var_int = catch_inputs$var,
                                    type = catch_inputs$type)

# NOTE: work out ATE procedure "by hand"
var_of_interest <- as.numeric(colData(methy_tmle_screened)[, catch_inputs$var])
if (class(var_of_interest) != "numeric") {
  var_of_interest <- as.numeric(var_of_interest)
}

# create clusters
methy_tmle_screened <- cluster_sites(methy_tmle = methy_tmle_screened)

# find all cases that have no missing values
cases_complete <- complete.cases(colData(methy_tmle_screened))

#...

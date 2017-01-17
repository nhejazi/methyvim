# R/`methadapt`

[![Travis-CI Build Status](https://travis-ci.org/nhejazi/methadapt.svg?branch=master)
](https://travis-ci.org/nhejazi/methadapt)
[![AppVeyor Build  Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/methadapt?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/methadapt)
[![Coverage Status](https://coveralls.io/repos/github/nhejazi/methadapt/badge.svg?branch=master)](https://coveralls.io/github/nhejazi/methadapt?branch=master)
[![RPKG](http://www.r-pkg.org/badges/version/methadapt)](http://www.r-pkg.org/pkg/methadapt)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)


`methadapt` is an R package that provides differential methylation analysis
and associated inference based on a statistical target parameter: the local
average treatment effect (L-ATE). The principal technique implemented here
relies on the use of targeted minimum loss-based estimation (TMLE) to analyze
and provide inference for a reduced set of genomic sites on which data is
typically made available by DNA methylation assays. In this procedure, multiple
testing corrections are made using a modified marginal Benjamini & Hochberg
step-up False Discovery Rate controlling procedure for multi-stage analyses
(FDR-MSA). While making allowance for the user to specify the reduced set of
genomic sites to be tested (e.g., those related to a particular gene of
interest), this package provides several screening procedures to facilitate
the initial reduction of genomic sites, including one based on the well-known
linear modeling method of the R package
[`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html) and
another based on data-adaptive statistical procedures for multiple testing,
implemented in the R package
[`data.adapt.multi.test`](https://github.com/wilsoncai1992/data.adapt.multi.test).

---

## Installation

- Install from GitHub: `devtools::install_github("nhejazi/methadapt", subdir = "pkg")`

---

## Principal References

* [Mark J. van der Laan & Sherri Rose. _Targeted Learning: Causal Inference for
    Observational and Experimental Data_,
    2011.](http://www.targetedlearningbook.com)

* [Catherine Tuglus & Mark J. van der Laan. "FDR controlling procedures for
    multi-stage analyses", _UC Berkeley Division of Biostatistics Working
    Paper Series_, 2008.](http://biostats.bepress.com/ucbbiostat/paper239/)

* [Weixin Cai, Nima S. Hejazi, & Alan E. Hubbard. "Data-adaptive statistics for
    multiple testing in high-dimensional settings", _In preparation_,
    2017.](https://www.overleaf.com/5660573pjjrxh#/25678897/)

* [Gordon K. Smyth. "Linear models and empirical Bayes methods for assessing
    differential expression in microarray experiments." _Statistical
    Applications in Genetics and Molecular Biology_, 3(1),
    2004.](http://www.statsci.org/smyth/pubs/ebayes.pdf)

---

## License

&copy; 2016-2017 [Nima S. Hejazi](http://nimahejazi.org), [Alan E.
Hubbard](http://sph.berkeley.edu/alan-hubbard), [Mark J. van der
Laan](https://www.stat.berkeley.edu/~laan/)

The contents of this repository are distributed under the MIT license. See file
`LICENSE` for details.

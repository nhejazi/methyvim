# R/`methyvim`

[![Travis-CI Build Status](https://travis-ci.org/nhejazi/methyvim.svg?branch=master)](https://travis-ci.org/nhejazi/methyvim)
[![AppVeyor Build  Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/methyvim?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/methyvim)
[![Coverage Status](https://img.shields.io/codecov/c/github/nhejazi/methyvim/master.svg)](https://codecov.io/github/nhejazi/methyvim?branch=master)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)


`methyvim` is an R package that provides facilities for differential methylation
analysis based on _variable importance measures_ (VIMs), a class of statistical
target parameters that arise in causal inference. The statistical techniques
implemented rely on the use of targeted minimum loss-based estimation (TMLE) to
assess two particular VIMs: (1) the "local" average treatment effect for
discrete exposures, and (2) a nonparametric variable importance measure for
continuous exposures. Within this framework, these methods allow differential
methylation effects to be estimated in an assumption-free manner, on a
pre-screened set of genomic sites measured by DNA methylation assays. As the
statistical algorithm uses a multi-stage approach, multiple testing corrections
are made using a modified marginal Benjamini & Hochberg step-up False Discovery
Rate controlling procedure for multi-stage analyses (FDR-MSA). In order to allow
the user significant flexibility with respect to the scientific questions posed,
the procedure estimates one of the appropriate VIMs for each of a reduced set of
genomic sites, using screening procedures to identify a reduced set from the
full data when a reduced set is not pre-specified. While making allowance for
the user to specify the reduced set of genomic sites to be tested (e.g., those
related to a particular gene of interest), this package comes equipped with
screening procedures to facilitate the initial reduction of genomic sites, the
two most noteworthy of these being a procedure based on the well-known R package
[`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html) and an
extension of data-adaptive statistical target parameters for multiple testing
(implemented in the R package
[`data.adapt.multi.test`](https://github.com/wilsoncai1992/data.adapt.multi.test)).

---

## Installation

- Install from GitHub: `devtools::install_github("nhejazi/methyvim")`

---

## Principal References

* [Mark J. van der Laan & Sherri Rose. _Targeted Learning: Causal Inference for
    Observational and Experimental Data_,
    2011.](http://www.targetedlearningbook.com)

* [Antoine Chambaz, Pierre Neuvial, & Mark J. van der Laan. "Estimation of a
    non-parametric variable importance measure of a continuous exposure",
    _Electronic Journal of Statistics_,
    2012.](http://www.math-info.univ-paris5.fr/~chambaz/Papiers/chambazNeuvialvanderLaan_EJS2012.pdf)

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

&copy; 2017 [Nima S. Hejazi](http://nimahejazi.org), [Alan E.
Hubbard](http://sph.berkeley.edu/alan-hubbard), [Mark J. van der
Laan](https://www.stat.berkeley.edu/~laan/)

The contents of this repository are distributed under the MIT license. See file
`LICENSE` for details.

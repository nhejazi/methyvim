
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`methyvim`
============

[![Travis-CI Build Status](https://travis-ci.org/nhejazi/methyvim.svg?branch=master)](https://travis-ci.org/nhejazi/methyvim) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/methyvim?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/methyvim) [![Coverage Status](https://img.shields.io/codecov/c/github/nhejazi/methyvim/master.svg)](https://codecov.io/github/nhejazi/methyvim?branch=master) [![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted learning of variable importance measures for the analysis of differential methylation experiments

**Author:** [Nima Hejazi](http://nimahejazi.org)

------------------------------------------------------------------------

Description
-----------

`methyvim` is an R package that provides facilities for differential methylation analysis based on *variable importance measures* (VIMs), a class of statistical target parameters that arise in causal inference.

The statistical methodology implemented computes targeted minimum loss-based estimates of well-studied variable importance measures:

1.  The average treatment effect (ATE) for discrete exposures/treatments,
2.  A nonparametric variable importance measure (NPVI) for continuous exposures (Chambaz, Neuvial, and van der Laan 2012).

These methods allow differential methylation effects to be quantified in a manner that is largely free of assumptions, especially of the variety exploited in parametric models. **The statistical algorithm consists in several major steps:**

1.  Pre-screening of genomic sites is used to isolate a subset of sites for which there is evidence of differential methylation. Targeted minimum loss-based estimates of VIMs are computed for this subset of sites only (for the sake of computational feasibility). Several screening approaches are available, adapting core routines from the following R packages: [`limma`](http://bioconductor.org/packages/release/bioc/html/limma.html), [`tmle.npvi`](https://CRAN.R-project.org/package=tmle.npvi), [`adaptest`](https://github.com/nhejazi/adaptest).
2.  Nonparametric VIMs are estimated for the specified parameter, currently using routines from the R packages [`tmle.npvi`](https://CRAN.R-project.org/package=tmle.npvi) to estimate the NPVI parameter, or [`tmle`](https://CRAN.R-project.org/package=tmle) to estimate the ATE.
3.  Since pre-screening is performed prior to estimating VIMs, we make use of a multiple testing correction uniquely suited to such settings, specifically we use the modified marginal Benjamini & Hochberg step-up False Discovery Rate controlling procedure for multi-stage analyses (FDR-MSA) (Tuglus and van der Laan 2009).

See van der Laan and Rose (2011) for a general description of the framework of targeted minimum loss-based estimation.

<!--
Note about shrinkage of influence curves, adapting @smyth2004linear.

More exposition here....
-->

------------------------------------------------------------------------

Installation
------------

<!--
For standard use, install from [Bioconductor](https://bioconductor.org):

```r
source("https://bioconductor.org/biocLite.R")
biocLite("methyvim")
```
-->
Install the most recent *stable release* from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("nhejazi/methyvim")
```

<!--
To contribute, install the _development version_ from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):


```r
devtools::install_github("nhejazi/methyvim", ref = "develop")
```
-->

------------------------------------------------------------------------

<!--
## Example

This is a basic example which shows you how to solve a common problem:


```r
## basic example code
```
-->
Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/nhejazi/methyvim/issues).

------------------------------------------------------------------------

Contributions
-------------

It is our hope that `methyvim` will grow to be widely adopted as a tool for the nonparametric assessment of variable importance in studies of differential methylation. To that end, contributions are very welcome, though we ask that interested contributors consult our [`contribution guidelines`](https://github.com/nhejazi/methyvim/blob/master/CONTRIBUTING.md) prior to submitting a pull request.

------------------------------------------------------------------------

Citation
--------

After using the `methyvim` R package, please cite both of the following:

        @article{hejazi2017methyvim,
          doi = {},
          url = {},
          year  = {2017},
          month = {},
          publisher={Faculty of 1000 Ltd},
          volume = {},
          author = {Hejazi, Nima S. and Hubbard, Alan E. and {van der Laan},
            Mark J.},
          title = {methyvim: Differential methylation analysis with targeted
            nonparametric variable importance measures},
          journal = {F1000Research}
        }

        @misc{hejazi2017methyvim-r,
          doi = {},
          url = {},
          howpublished = {\url{https://github.com/nhejazi/methyvim}},
          year  = {2017},
          author = {Hejazi, Nima S. and Hubbard, Alan E. and {van der Laan},
            Mark J.},
          title = {methyvim: Differential methylation analysis with targeted
            minimum loss-based estimates of variable importance measures}
        }

------------------------------------------------------------------------

Related
-------

-   [R/`methyvimData`](https://github.com/nhejazi/methyvimData) - R package with sample experimental DNA methylation data for use as an example with this analysis package.

------------------------------------------------------------------------

License
-------

© 2017 [Nima S. Hejazi](http://nimahejazi.org), [Alan E. Hubbard](http://sph.berkeley.edu/alan-hubbard), [Mark J. van der Laan](https://www.stat.berkeley.edu/~laan/)

The contents of this repository are distributed under the MIT license. See file `LICENSE` for details.

------------------------------------------------------------------------

References
----------

Chambaz, Antoine, Pierre Neuvial, and Mark J van der Laan. 2012. “Estimation of a Non-Parametric Variable Importance Measure of a Continuous Exposure.” *Electronic Journal of Statistics* 6. NIH Public Access: 1059.

Tuglus, Catherine, and Mark J. van der Laan. 2009. “Modified FDR Controlling Procedure for Multi-Stage Analyses.” *Statistical Applications in Genetics and Molecular Biology* 8 (1). Walter de Gruyter: 1–15. doi:[10.2202/1544-6115.1397](https://doi.org/10.2202/1544-6115.1397).

van der Laan, Mark J., and Sherri Rose. 2011. *Targeted Learning: Causal Inference for Observational and Experimental Data*. Springer Science & Business Media.

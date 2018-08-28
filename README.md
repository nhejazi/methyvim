
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`methyvim`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/methyvim.svg?branch=master)](https://travis-ci.org/nhejazi/methyvim)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/methyvim?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/methyvim)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/methyvim/master.svg)](https://codecov.io/github/nhejazi/methyvim?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/methyvim.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/methyvim)
[![Bioc
Time](http://bioconductor.org/shields/years-in-bioc/methyvim.svg)](https://bioconductor.org/packages/release/bioc/html/methyvim.html)
[![Bioc
Downloads](http://bioconductor.org/shields/downloads/methyvim.svg)](https://bioconductor.org/packages/release/bioc/html/methyvim.html)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1401298.svg)](https://doi.org/10.5281/zenodo.1401298)

> Targeted, Robust, and Model-free Differential Methylation Analysis

**Author:** [Nima Hejazi](https://nimahejazi.org)

-----

## What’s `methyvim`?

`methyvim` is an R package that provides facilities for differential
methylation analysis based on *variable importance measures* (VIMs), a
class of statistically estimable target parameters that arise in causal
inference.

The statistical methodology implemented computes targeted minimum
loss-based estimates of several well-characterized variable importance
measures:

For discrete-valued treatments or exposures:

  - The *average treatment effect* (ATE): The effect of a binary
    exposure or treatment on the observed methylation at a target CpG
    site is estimated, controlling for the observed methylation at all
    other CpG sites in the same neighborhood as the target site, based
    on an additive form. In particular, the parameter estimate
    represents the **additive difference** in methylation that would
    have been observed at the target site had all observations received
    the treatment versus the scenario in which none received the
    treatment.

  - The *relative risk* (RR): The effect of a binary exposure or
    treatment on the observed methylation at a target CpG site is
    estimated, controlling for the observed methylation at all other CpG
    sites in the same neighborhood as the target site, based on an
    geometric form. In particular, the parameter estimate represents the
    **multiplicative difference** in methylation that would have been
    observed at the target site had all observations received the
    treatment versus the scenario in which none received the treatment.

For continous-valued treatments or exposures:

  - A *nonparametric variable importance measure* (NPVI) (Chambaz,
    Neuvial, and van der Laan 2012): The effect of continous-valued
    exposure or treatment (the observed methylation at a target CpG
    site) on an outcome of interest is estimated, controlling for the
    observed methylation at all other CpG sites in the same neighborhood
    as the target (treatment) site, based on a parameter that compares
    values of the treatment against a reference value taken to be the
    null. In particular, the implementation provided is designed to
    assess the effect of differential methylation at the target CpG site
    on a (typically) phenotype-level outcome of interest (e.g.,
    survival), in effect providing an nonparametric evaluation of the
    impact of methylation at the target site on said outcome.

*In all cases, an estimator of the target parameter is constructed via
targeted minimum loss-based estimation.*

These methods allow differential methylation effects to be quantified in
a manner that is largely free of assumptions, especially of the variety
exploited in parametric models. **The statistical algorithm consists in
several major steps:**

1.  Pre-screening of genomic sites is used to isolate a subset of sites
    for which there is cursory evidence of differential methylation. For
    the sake of computational feasibility, targeted minimum loss-based
    estimates of VIMs are computed only for this subset of sites.
    Several screening approaches are available, adapting core routines
    from the following R packages:
    [`limma`](http://bioconductor.org/packages/release/bioc/html/limma.html),
    [`tmle.npvi`](https://CRAN.R-project.org/package=tmle.npvi).
2.  Nonparametric VIMs are estimated for the specified parameter,
    currently adapting routines from the
    [`tmle.npvi`](https://CRAN.R-project.org/package=tmle.npvi) and
    [`tmle`](https://CRAN.R-project.org/package=tmle) R packages.
3.  Since pre-screening is performed prior to estimating VIMs, we make
    use of a multiple testing correction uniquely suited to such
    settings. Due to the multiple testing nature of the estimation
    problem, a variant of the Benjamini & Hochberg procedure for
    controlling the False Discovery Rate (FDR) is applied (Benjamini and
    Hochberg 1995). Specifically, we apply the modified marginal
    Benjamini & Hochberg step-up False Discovery Rate controlling
    procedure for multi-stage analyses (FDR-MSA) (Tuglus and van der
    Laan 2009).

For a general discussion of the framework of targeted minimum loss-based
estimation, the many applications of this methodology, and the role the
framework plays in statistical causal inference, the recommended
references are van der Laan and Rose (2011) and van der Laan and Rose
(2018). Hernan and Robins (2018) and Pearl (2009) may be of interest to
those desiring a more general introduction to statistical causal
inference.

-----

## Installation

For standard use, install from
[Bioconductor](https://bioconductor.org/packages/methyvim) using
[`BiocManager`](https://CRAN.R-project.org/package=BiocManager):

``` r
if (!("BiocManager" %in% installed.packages())) {
  install.packages("BiocManager")
}
BiocManager::install("methyvim")
```

To contribute, install the bleeding-edge *development version* from
GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("nhejazi/methyvim")
```

Current and prior [Bioconductor](https://bioconductor.org) releases are
available under branches with numbers prefixed by “RELEASE\_”. For
example, to install the version of this package available via
Bioconductor 3.6, use

``` r
devtools::install_github("nhejazi/methyvim", ref = "RELEASE_3_6")
```

-----

## Example

For details on how to best use the `methyvim` R package, please consult
the most recent [package
vignette](http://bioconductor.org/packages/release/bioc/vignettes/methyvim/inst/doc/using_methyvim.html)
available through the [Bioconductor
project](https://bioconductor.org/packages/methvim).

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/methyvim/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/methyvim/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `methyvim` R package, please cite the following:

``` 
    @article{hejazi2018methyvim,
      doi = {},
      url = {},
      year  = {2018},
      month = {},
      publisher = {},
      volume = {},
      author = {Hejazi, Nima S and Phillips, Rachael V and Hubbard, Alan E
        and {van der Laan}, Mark J},
      title = {methyvim: Targeted, Robust, and Model-free Differential
        Methylation Analysis in {R}},
      journal = {}
    }
```

-----

## Related

  - [R/`methyvimData`](https://github.com/nhejazi/methyvimData) - R
    package with sample experimental DNA methylation data for use as an
    example with this analysis package.

-----

## Funding

The development of this software was supported in part through grants
from the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=),
[R01
ES021369-05](https://projectreporter.nih.gov/project_info_description.cfm?aid=9210551&icde=37849782&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=),
and [P42
ES004705-29](https://projectreporter.nih.gov/project_info_details.cfm?aid=9260357&map=y).

-----

## License

© 2017-2018 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See file `LICENSE` for details.

-----

## References

<div id="refs" class="references">

<div id="ref-benjamini1995controlling">

Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False
Discovery Rate: A Practical and Powerful Approach to Multiple Testing.”
*Journal of the Royal Statistical Society. Series B (Methodological)*.
JSTOR, 289–300.

</div>

<div id="ref-chambaz2012estimation">

Chambaz, Antoine, Pierre Neuvial, and Mark J van der Laan. 2012.
“Estimation of a Non-Parametric Variable Importance Measure of a
Continuous Exposure.” *Electronic Journal of Statistics* 6. NIH Public
Access: 1059.

</div>

<div id="ref-hernan2018causal">

Hernan, Miguel A, and James M Robins. 2018. *Causal Inference*. Chapman
& Hall / CRC Texts in Statistical Science. Taylor & Francis.

</div>

<div id="ref-pearl2009causality">

Pearl, Judea. 2009. *Causality: Models, Reasoning, and Inference*.
Cambridge University Press.

</div>

<div id="ref-tuglus2009modified">

Tuglus, Catherine, and Mark J van der Laan. 2009. “Modified FDR
Controlling Procedure for Multi-Stage Analyses.” *Statistical
Applications in Genetics and Molecular Biology* 8 (1). Walter de
Gruyter: 1–15. <https://doi.org/10.2202/1544-6115.1397>.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vdl2018targeted">

———. 2018. *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*. Springer Science & Business Media.

</div>

</div>

---
title: '`methadapt`: differential methylation analysis and local inference with
        targeted minimum loss-based estimation'
tags:
  - Targeted Minimum Loss-Based Estimation
  - multiple testing
  - causal inference
  - DNA methylation
  - bioinformatics
  - genomics
authors:
 - name: Nima S. Hejazi
   orcid: 0000-0002-7127-2789
   affiliation: 1
 - name: Alan E. Hubbard
   orcid: 0000-0002-3769-0127
   affiliation: 1
 - name: Mark J. van der Laan
   orcid:
   affiliation: 1
affiliations:
 - name: Division of Biostatistics, University of California, Berkeley
   index: 1
date: 01 May 2017
bibliography: paper.bib
---

# Summary

`methadapt` is an R package that provides differential methylation analysis
and associated inference based on a statistical target parameter: the local
average treatment effect (L-ATE). The principal technique implemented here
relies on the use of targeted minimum loss-based estimation (TMLE) to analyze
and provide inference for a reduced set of genomic sites on which data is
typically made available by DNA methylation assays. In this procedure, multiple
testing corrections are made using a modified marginal Benjamini & Hochberg
step-up False Discovery Rate controlling procedure for multi- stage analyses
(FDR-MSA), proposed by @tuglus2008fdr. Although it is possible for the user to
specify the reduced set of genomic sites to be tested (e.g., those related to a
particular gene of interest), this package provides several screening
procedures to facilitate the initial reduction of genomic sites as well,
including one based on the well-known linear modeling method of the R package
LIMMA (@smyth2005limma) and another based on data-adaptive statistical
procedures for multiple testing (@cai207data).

\newpage

# References

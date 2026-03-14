# oligodendrogenesis-Parkinson-s-disease
Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease
## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)

# Overview

This repository provides the computational framework and analytical pipelines used in our study: “Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease.” By integrating multi-omics analyses with cell type-specific genetic manipulation to identify OLs as critical drivers of PD pathogenesis. Partitioned heritability analysis revealed significant enrichment of PD genetic risk in OLs, and single-nucleus RNA sequencing demonstrated consistent OL expansion in the SN of human postmortem tissue and two independent mouse models. Transcriptomic integration identified OLs as the primary cellular source of ROS in the PD SN.

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation, and usage of the `lolR` package on many real and simulated data examples.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.


# System Requirements

Bulk RNA-seq analysis: R version 4.4.1, Bioconductor version 3.19, DESeq2 version 1.44.0, enrichplot version 1.24.4, clusterProfiler version 4.12.6, msigdbr version 24.1.0, fgsea version 1.30.0, limma version 3.60.6.
ATAC-seq analysis: R version 4.5.1, Signac version 1.14.9002.
ScRNA-seq analysis: R version 4.5.1, SeuratObject version 5.3.0, harmony version 1.2.3, samtools version 1.21, cellranger version 8.0.0, FastQC version 0.12.1, hisat2 version 2.2.1.
Partitioned Heritability Analysis: LDSC version 1.0.1 (https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format ).


# Instructions for Use


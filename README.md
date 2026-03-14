# oligodendrogenesis-Parkinson-s-disease
Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease
## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Instructions for Use](#Instructions-for-Use)
- [Citation Database](#Citation-Database)

# Overview

This repository provides the computational framework and analytical pipelines used in our study: “Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease.” By integrating multi-omics analyses with cell type-specific genetic manipulation to identify OLs as critical drivers of PD pathogenesis. Partitioned heritability analysis revealed significant enrichment of PD genetic risk in OLs, and single-nucleus RNA sequencing demonstrated consistent OL expansion in the SN of human postmortem tissue and two independent mouse models. Transcriptomic integration identified OLs as the primary cellular source of ROS in the PD SN.

# Repo Contents

- [LDSC](./figure1_a&b_ldsc): `LDSC` code.
- [EWCE](./figure1_c_EWCE): `EWCE` code.
- [scRNA-seq](./figure2&figure3_i-l_scRNA-seq_analysis): `scRNA-seq` code.
- [Bulk RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq): `Bulk RNA-seq` code.
- [AUCell](./figure3_d&h_AUCell): `AUCell` code.


# System Requirements

Bulk RNA-seq analysis: R version 4.4.1, Bioconductor version 3.19, DESeq2 version 1.44.0, enrichplot version 1.24.4, clusterProfiler version 4.12.6, msigdbr version 24.1.0, fgsea version 1.30.0, limma version 3.60.6.
ATAC-seq analysis: R version 4.5.1, Signac version 1.14.9002.
ScRNA-seq analysis: R version 4.5.1, SeuratObject version 5.3.0, harmony version 1.2.3, samtools version 1.21, cellranger version 8.0.0, FastQC version 0.12.1, hisat2 version 2.2.1.
Partitioned Heritability Analysis: LDSC version 1.0.1 (https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format ).


# Instructions for Use
*LDSC在此文档文件中能够正常运行，如有报错，请检查配置环境。
- [LDSC](./figure1_a&b_ldsc): 
  本文件包含LDSC运行的过程代码，包含从单细胞数据中提取差异表达基因、基线数据生成与遗传力计算和p值计算的一系列代码。总共有三个步骤，步骤1和步骤3为在`R`环境下运行的代码。步骤2为在`Linux`环境下运行的LDSC软件包程序。
- [EWCE](./figure1_c_EWCE):
  本文件包含对单细胞基因集的打分代码，代码中包含数据预处理与特异性计算，EWCE分析函数构建，差异基因提取与基于背景基因集的打分。
- [scRNA-seq](./figure2&figure3_i-l_scRNA-seq_analysis):
  本文件包含对五组单细胞数据的处理代码，整合数据的前期处理单独放置于[整合数据前处理](./figure2&figure3_i-l_scRNA-seq_analysis/数据整合前处理)中。CNP0000892的单细胞数据的处理过程与MPTP数据的处理过程及整合数据的处理过程分别命名为[CNP0000892](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_d&e&f(CNP).R)，[MPTP](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_g&h&i(MPTP).R)，[整合数据](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_a&b&c(整合数据).R)。整合数据亚聚类代码与亚聚类的GO通路富集分析流程保存在同目录下的figure3_i&j(亚聚类).R与figure3_k&l(GO).R中
- [Bulk_RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq):
  本文件包含对两组简单RNA测序数据的处理代码，差异基因分析流程文件放置于[前处理](./figure3_a-c&e-g_Bulk_RNA-seq/前处理)中。火山图、热图与富集通路分析的GO图的运行流程代码均放置于[Bulk_RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq)目录下。
- [AUCell](./figure3_d&h_AUCell):
  本文件包含对两组简单RNA测序数据在整合后的单细胞数据的打分代码。
  
# Citation Database

Four PD GWAS summary statistics datasets derived primarily from individuals of European ancestry were compiled. These included: Leal et al., 2025 (doi: 10.1101/2025.07.18.25331793) ; GP2 Parkinson’s Disease GWAS, European ancestry, 2025; Kim et al., 2024 (doi: 10.1038/s41588-023-01584-8);  and Nalls et al., 2019 (doi: 10.1016/S1474-4422(19)30320-5).
6 human PD datasets (GSE243639, GSE178265, GSE193688, GSE147672, GSE7621 and GSE49036) were obtained from the GEO database. The MPTP mouse dataset was retrieved from Figshare (identifier: https://figshare.com/articles/dataset/Single cell_sequencing_screening_of_biomarkers_associated_with_Lipophagy_in_the_MPTP_model_mice_of_Parkinson_s_disease/26724526/1?file=48569005), and the α-synuclein A53T transgenic mouse dataset was obtained from the CNGBdb (CNSA; accession code CNP0000892, https://db.cngb.org/cnsa/).

# Bulk RNA-seq Differential Expression and Visualization

This repository contains an R-based workflow for bulk RNA-seq analysis focused on comparing control and treated samples. 
In addition, the workflow includes downstream visualization and exploratory analysis such as volcano plots and fold-change correlation.

## Overview

This workflow is designed to:

- prepare count-based RNA-seq data for analysis
- define control and treated sample groups
- perform differential expression analysis using `edgeR`
- identify significantly upregulated and downregulated genes
- generate volcano plots for visualization
- explore correlation across log fold-change comparisons

## Workflow Summary

### 1. Data Preparation
The script:
- reads raw STAR/RSEM count data
- sets Ensembl IDs as row names
- filters the dataset to retain relevant samples
- defines comparison groups for downstream analysis

### 2. Differential Expression Analysis
Differential expression is performed using `edgeR`, including:
- creation of a `DGEList` object
- filtering lowly expressed genes with `filterByExpr()`
- normalization with `calcNormFactors()`
- dispersion estimation
- statistical testing with `exactTest()`
- false discovery rate (FDR) correction

### 3. Significant Gene Filtering
The workflow extracts significantly differentially expressed genes using thresholds based on:
- adjusted p-value (FDR)
- log2 fold-change cutoffs

### 4. Volcano Plot Visualization
The script generates a volcano plot to highlight:
- upregulated genes
- downregulated genes
- non-significant genes

Gene symbols are mapped from Ensembl IDs.

### 5. Correlation Analysis
The workflow also includes correlation analysis of log fold-change values across comparisons to assess agreement between expression changes under different conditions.

## Tools and Packages

- R
- edgeR
- ggplot2
- ggrepel

## Input Requirements

This workflow assumes:
- count-based bulk RNA-seq input data
- samples organized in a count matrix format
- Ensembl gene identifiers
- sample names containing enough information to define groups

## Output

The workflow can generate:
- differential expression result tables
- filtered significant gene tables
- volcano plot figures
- correlation plots and summary statistics

## Notes

- Some file paths in the original script may be project-specific and should be updated before reuse.
- This workflow reflects a project-specific implementation and can be adapted into a more general reusable pipeline.
- Gene annotation is added by mapping Ensembl IDs to gene symbols through an external reference file.

## Repository Purpose

This repository demonstrates practical experience in:
- bulk RNA-seq analysis
- differential expression testing
- transcriptomic data visualization
- reproducible analysis development in R

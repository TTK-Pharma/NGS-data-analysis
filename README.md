![R](https://img.shields.io/badge/R-Programming-blue)
![Seurat](https://img.shields.io/badge/Seurat-scRNAseq-green)
![Field](https://img.shields.io/badge/Field-Bioinformatics-orange)

---
#Project preview
UMAP visualization of single-cell RNA sequencing data showing distinct cellular clusters identified from a glioblastoma dataset using the Seurat analysis workflow.
![Annotated UMAP](plots/Annotated-UMAp.png)

## Project overview
This project performs **single-cell RNA sequencing (scRNA-seq) analysis** of a glioblastoma dataset obtained from **10x Genomics**.

The analysis was conducted using the **Seurat package in R** to investigate cellular heterogeneity within the tumor microenvironment. The workflow includes preprocessing, quality control, normalization, dimensionality reduction, clustering, marker gene identification, and cluster annotation.

Single-cell RNA sequencing enables the identification of transcriptionally distinct cell populations within tumors, providing insights into tumor biology and immune interactions.

--
## Dataset
Dataset used:
**Parent_SC3v3_Human_Glioblastoma_raw_feature_bc_matrix.h5**
Source: 
10x Genomics public datasets
Technology:
Single-cell RNA Sequencing (scRNA-seq)
---

## Analysis Pipeline

The scRNA-seq analysis workflow includes the following steps:
Raw 10x Genomics Data (.h5)
↓
Quality Control Filtering
↓
Mitochondrial Gene Percentage Calculation
↓
Normalization
↓
Highly Variable Gene Identification
↓
Scaling of Gene Expression
↓
Principal Component Analysis (PCA)
↓
Cell Clustering
↓
UMAP Visualization
↓
Marker Gene Identification
↓
Cluster Annotation

---

## Tools and Packages Used

The analysis was performed using the following R packages:

- Seurat
- dplyr
- ggplot2
- patchwork

---

## Repository Structure

plots/
├── raw_plot.png
├── variableplot.png
├── elbowplot.png
├── dimplot.png
├── expressed_gene5.png
└── Annotated-UMAP.png

mainfile.R
markersByclusters.csv
README.md

---

## Reproducibility

To reproduce the analysis:

1. Install the required R packages
2. Download the dataset from 10x Genomics
3. Run the analysis script

---

## Environment

R version: 4.5.2  
Operating System: Windows 11

---

## Author

**Thirukumaran**

Bioinformatics | Computational Biology | Single-Cell RNA-seq Analysis

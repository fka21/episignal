# Differential Gene Expression Analysis with limma-voom and Wilcoxon Rank-Sum Test
 
This README provides an overview of the workflow for differential gene expression analysis on RNA-seq data, comparing normal and apical-out polarities in pancreatic organoid samples. The analysis combines statistical and enrichment approaches, with outputs for visualization saved to the `output/` directory.

## Overview
This repository contains scripts which conduct differential gene expression analysis on RNA-seq data using both **limma-voom** and **Wilcoxon Rank-Sum Test** methods. It visualizes the results using boxplots, PCA, heatmaps, and pathway enrichment plots. The gene IDs are in Ensembl format, and pathway analysis is performed using GO and KEGG human databases.

### Key Features
1. **Differential Expression Analysis**: 
   - *limma-voom*: Empirical Bayes moderated t-tests for robust gene expression analysis.
   - *Wilcoxon Rank-Sum Test*: A non-parametric paired test suitable for small sample sizes.

2. **Pathway Enrichment Analysis**: 
   - GO and KEGG pathway analysis, conducted with **clusterProfiler** and **org.Hs.eg.db**.

3. **Visualization and Output**: 
   - Various plots, including PCA, heatmaps, and enrichment dot plots, provide insights into the biological significance of the results.
   - Visualizations are stored in the `output/` directory for easy reference.

4. **Machine Learning (PyCaret)**: 
   - Machine learning approaches, were also explored to assess data variability and feature importance.

## Dependencies
The script is designed to install missing dependencies and utilizes a path-insensitive approach with relative paths. Required packages include:
- `tidyverse`, `ggfortify`, `limma`, `here`, `VennDiagram`, `clusterProfiler`, `org.Hs.eg.db`, `pheatmap`

## Data and Setup
- **Data**: RNA-seq TPM values are used for differential expression analysis.
- **Sample Categories**:
  - Normal polarity samples: `CG.in_S31`, `Nav.in_S36`, `QB.in_S33`
  - Apical-out organoids: `CG.out_S32`, `Nav.out_S1`, `QB.out_S34`
- **Directory Structure**: Creates an `output/` directory for saving results.

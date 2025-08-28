RNA-seq Analysis from Raw Counts

This repository provides a complete pipeline for RNA-seq analysis, from raw count data to differential expression and visualization.

1. Input
raw_counts.csv – Gene expression matrix (rows = genes, columns = samples) from STAR/HISAT2 + featureCounts/HTSeq.

MetaData.csv – Sample annotation table (condition, tissue type, replicate).

2. Analysis
RNAseqAnalysisCode.R implements:

Data import, QC, and DESeq2 normalization

PCA and clustering

Differential expression (DEG) analysis across conditions

Visualization: volcano plots, MA plots, heatmaps



3. Outputs

DeSeqResults.csv

Upregulated_* / Downregulated_* CSVs per comparison

Normalized_Counts.csv – DESeq2-normalized matrix


PCA_data.csv – PCA coordinates


Plots (Output/DESeq2_Plots/)

Heatmap (top 20 DEGs), MA plots, PCA, Volcano plots

4. Usage
Install Requirements

install.packages(c("tidyverse", "pheatmap", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ComplexHeatmap"))

 raw_counts.csv and MetaData.csv are in Input_data/


Update file paths if needed

Execute script → Results in Output/

5. Interpretation
Upregulated: Higher expression in first condition

Downregulated: Lower expression in first condition

Volcano plot: log₂FC vs. –log₁₀(p)

Heatmap: Clustering of top DEGs

Citation:
Dawadi P., Roy M., Pokharel B., Shrestha A., Niraula N., Naeem A., Miura S., Nepal S*. Department of Biology, University of Mississippi, MS 38677, USA. Under review. 




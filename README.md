
# RNA-seq Differential Expression (DESeq2) — Code Description

This document explains each step of your R script for RNA-seq differential expression analysis, what it expects as input, what it produces as output.

---

## 1) Libraries & Setup

**What it does**
- Loads analysis and plotting libraries: `DESeq2`, `ggplot2`, `pheatmap`, `RColorBrewer`, `ashr`, `readr`, `dplyr`, `matrixStats`, `vegan`.
- Sets the working directory and creates an output folder for plots.

**Key lines**
```r
library(DESeq2); library(ggplot2); library(pheatmap); library(RColorBrewer)
library(ashr); library(readr); library(dplyr); library(matrixStats); library(vegan)

setwd("C:/Users/newfaculty/Desktop/Bioinformatics_Club")
plot_dir <- "DESeq2_Plots"
dir.create(plot_dir, showWarnings = FALSE)
```

**Inputs**: None yet.  
**Outputs**: Creates `DESeq2_Plots/` folder.  
**Tweak**: Change `setwd(...)` and `plot_dir` as needed.

---

## 2) Load Inputs & Align Samples

**What it does**
- Reads count matrix and metadata.
- Extracts the grouping column `biopsy_site` from metadata.
- Ensures column order of counts matches the row order of metadata (critical for DESeq2).

**Key lines**
```r
counts <- read.csv("raw_counts.csv", row.names = 1)
meta   <- read.csv("MetaData.csv", row.names = 1)

colnames(meta) <- trimws(colnames(meta))
meta$biopsy_site <- meta$Factor.Value.biopsy.site.

counts <- counts[, rownames(meta)]
stopifnot(all(colnames(counts) == rownames(meta)))
```

**Inputs**: `raw_counts.csv` (genes x samples), `MetaData.csv` (samples x attributes).  
**Outputs**: In-memory objects `counts`, `meta`.  
**Tweak**: If the biopsy-site column has a different name, adjust `meta$Factor.Value.biopsy.site.`.

---

## 3) Filter Low-Count Genes

**What it does**
- Drops genes with fewer than 10 counts in at least 2 samples (noise reduction).

**Key lines**
```r
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]
message("Genes retained after filtering: ", nrow(counts))
```

**Inputs**: `counts`.  
**Outputs**: Filtered `counts`.  
**Tweak**: Thresholds (10 counts, ≥2 samples) can be changed depending on library size and design.

---

## 4) Define Factor Levels (Reference Group)

**What it does**
- Orders levels so that `normal` is the reference group for contrasts.
- Levels used: `normal`, `primary tumor`, `colorectal cancer metastatic in the liver`.

**Key lines**
```r
meta$biopsy_site <- factor(
  meta$biopsy_site,
  levels = c("normal", "primary tumor", "colorectal cancer metastatic in the liver")
)
```

**Tweak**: Ensure your metadata values exactly match these strings or edit accordingly.

---

## 5) Build DESeq2 Dataset & Run DESeq

**What it does**
- Constructs the `DESeqDataSet` with design `~ biopsy_site` and runs dispersion estimation and model fitting.

**Key lines**
```r
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ biopsy_site)

dds <- DESeq(dds)
```

**Outputs**: `dds` object with fitted model.  
**Note**: Counts are rounded to integers for DESeq2.

---

## 6) Variance-Stabilizing Transform (VST)

**What it does**
- Produces normalized, variance-stabilized expression for visualization and ordination.
- Saves normalized matrix to CSV.

**Key lines**
```r
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
write.csv(norm_counts, file = "Normalized_Counts.csv")
```

**Outputs**: `Normalized_Counts.csv` (VST values).  
**Tweak**: `blind = FALSE` respects design; set `TRUE` for purely exploratory, design-agnostic transform.

---

## 7) Global Results Table (Default Contrast)

**What it does**
- Extracts default results (last vs first level of `biopsy_site`), orders by FDR, saves to CSV.

**Key lines**
```r
res <- results(dds)
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = "DESeq2_results.csv")
```

**Outputs**: `DESeq2_results.csv` (all genes, default contrast).

---

## 8) PCA Plot of Samples

**What it does**
- Computes PCA from `vsd`, annotates by `biopsy_site`, and saves a publication-style figure.

**Key lines**
```r
pcaData <- plotPCA(vsd, intgroup = "biopsy_site", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(file.path(plot_dir, "PCA_plot.png"), width = 7, height = 6, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, color = biopsy_site)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples by Biopsy Site") +
  theme_minimal(base_size = 20) +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
dev.off()

write.csv(pcaData, file = "PCA_data.csv", row.names = TRUE)
```

**Outputs**: `DESeq2_Plots/PCA_plot.png`, `PCA_data.csv`.

---

## 9) Differential Expression: Primary vs Normal

**What it does**
- Computes contrast primary tumor vs normal.
- Applies `ashr::lfcShrink` for stabilized log2FC.
- Saves ordered table and MA/Volcano plots.

**Key lines**
```r
res_primary <- results(dds, contrast = c("biopsy_site", "primary tumor", "normal"))
res_primary <- lfcShrink(dds, coef = "biopsy_site_primary.tumor_vs_normal", type = "ashr", res = res_primary)
res_primary_ordered <- res_primary[order(res_primary$pvalue), ]
write.csv(as.data.frame(res_primary_ordered), "Primary_vs_Normal_DEGs.csv")

# MA plot
png(file.path(plot_dir, "MA_plot_Primary_vs_Normal.png"), width = 15, height = 14, units = "in", res = 300)
par(cex = 1.5)
plotMA(res_primary, ylim = c(-5, 5), main = "MA Plot: Primary Tumor vs Normal")
dev.off()

# Volcano plot
df_primary <- as.data.frame(res_primary)
df_primary <- df_primary[is.finite(df_primary$log2FoldChange) & is.finite(df_primary$padj), ]
df_primary$threshold <- ifelse(df_primary$padj < 0.05 & df_primary$log2FoldChange > 1, "Upregulated",
                        ifelse(df_primary$padj < 0.05 & df_primary$log2FoldChange < -1, "Downregulated", "NS"))
# (plot code follows…)
```

**Outputs**
- Table: `Primary_vs_Normal_DEGs.csv`
- Plots: `DESeq2_Plots/MA_plot_Primary_vs_Normal.png`, `DESeq2_Plots/Volcano_Primary_vs_Normal.png`

**Tweak**
- Volcano thresholds: `padj < 0.05` and `|log2FC| > 1`. Adjust for your study (e.g., stricter FC or FDR).

---

## 10) Differential Expression: Liver Mets vs Normal

**What it does**
- Computes contrast liver metastasis vs normal with shrinkage, saves table and MA/Volcano plots.

**Key lines**
```r
res_liver <- results(dds, contrast = c("biopsy_site", "colorectal cancer metastatic in the liver", "normal"))
res_liver <- lfcShrink(dds, coef = "biopsy_site_colorectal.cancer.metastatic.in.the.liver_vs_normal", type = "ashr", res = res_liver)
res_liver_ordered <- res_liver[order(res_liver$pvalue), ]
write.csv(as.data.frame(res_liver_ordered), "LiverMets_vs_Normal_DEGs.csv")
# (MA and Volcano plot code follows…)
```

**Outputs**
- Table: `LiverMets_vs_Normal_DEGs.csv`
- Plots: `DESeq2_Plots/MA_plot_Liver_vs_Normal.png`, `DESeq2_Plots/Volcano_Liver_vs_Normal.png`

---

## 11) Heatmaps: Top 20 Variable Genes

**What it does**
- Selects top 20 most variable genes (via row variance on VST matrix).
- Produces **two heatmaps**:
  1) **Unscaled (VST values)** using a sequential palette with winsorized color breaks (2nd–98th percentiles).
  2) **Row Z-score scaled** heatmap using a diverging palette; caps z-scores at ±2.5 for robust contrast.

**Key lines**
```r
n_top <- 20; z_cap <- 2.5
mat_all <- assay(vsd)
sel_idx <- head(order(rowVars(mat_all), decreasing = TRUE), n_top)
heatmap_data <- as.matrix(mat_all[sel_idx, , drop = FALSE])

annotation_col <- data.frame(Biopsy = meta$biopsy_site)
rownames(annotation_col) <- colnames(heatmap_data)

# Unscaled with winsorized breaks
q <- quantile(heatmap_data, probs = c(0.02, 0.98), na.rm = TRUE)
breaks_cont <- seq(q[1], q[2], length.out = 102)

# Z-score per row with capping
z_mat <- t(scale(t(heatmap_data)))
z_mat[!is.finite(z_mat)] <- 0
z_mat[z_mat > z_cap] <- z_cap; z_mat[z_mat < -z_cap] <- -z_cap
breaks_div <- seq(-z_cap, z_cap, length.out = 102)
```

**Outputs**
- `DESeq2_Plots/Heatmap_Top20Genes_unscaled.png`
- `DESeq2_Plots/Heatmap_Top20Genes_Zscore.png`

**Tweak**
- Change `n_top`, `z_cap`, clustering distance/method, and palettes to suit your figure style.

---

## 12) Export DEG Subsets

**What it does**
- Saves significantly up- and downregulated genes for each contrast (padj < 0.05 and |log2FC| > 1).

**Key lines**
```r
up_primary   <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange >  1)
down_primary <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_primary,   "Upregulated_Primary_vs_Normal.csv")
write.csv(down_primary, "Downregulated_Primary_vs_Normal.csv")

up_liver   <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange >  1)
down_liver <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_liver,   "Upregulated_Liver_vs_Normal.csv")
write.csv(down_liver, "Downregulated_Liver_vs_Normal.csv")
```

**Outputs**
- `Upregulated_Primary_vs_Normal.csv`
- `Downregulated_Primary_vs_Normal.csv`
- `Upregulated_Liver_vs_Normal.csv`
- `Downregulated_Liver_vs_Normal.csv`

---

## 13) File Output Summary (Quick Checklist)

- **Tables**
  - `DESeq2_results.csv`
  - `Normalized_Counts.csv`
  - `PCA_data.csv`
  - `Primary_vs_Normal_DEGs.csv`
  - `LiverMets_vs_Normal_DEGs.csv`
  - `Upregulated_Primary_vs_Normal.csv`
  - `Downregulated_Primary_vs_Normal.csv`
  - `Upregulated_Liver_vs_Normal.csv`
  - `Downregulated_Liver_vs_Normal.csv`

- **Figures (folder: `DESeq2_Plots/`)**
  - `PCA_plot.png`
  - `MA_plot_Primary_vs_Normal.png`
  - `Volcano_Primary_vs_Normal.png`
  - `MA_plot_Liver_vs_Normal.png`
  - `Volcano_Liver_vs_Normal.png`
  - `Heatmap_Top20Genes_unscaled.png`
  - `Heatmap_Top20Genes_Zscore.png`

---


## 16) References

Dawadi P., Pokharel B., Shrestha A., Niraula D., Naeema A., Miura S., Roy M.*, Nepal S.* (2025). From bench to bytes: A practical guide to RNA sequencing data analysis. Under review.   



# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ashr)
library(readr)
library(dplyr)
library(matrixStats)
library(vegan)


# Set working directory and output directory
setwd("C:/Users/newfaculty/Desktop/Bioinformatics_Club")
plot_dir <- "DESeq2_Plots"
dir.create(plot_dir, showWarnings = FALSE)

# Read input files
counts <- read.csv("raw_counts.csv", row.names = 1)
meta   <- read.csv("MetaData.csv", row.names = 1)

# Clean column names and extract condition info
colnames(meta) <- trimws(colnames(meta))
meta$biopsy_site <- meta$Factor.Value.biopsy.site.

# Ensure sample order matches
counts <- counts[, rownames(meta)]
stopifnot(all(colnames(counts) == rownames(meta)))

# Filter low-count genes: keep genes with ≥10 counts in ≥2 samples
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]
message("Genes retained after filtering: ", nrow(counts))

# Set reference level for biopsy_site
meta$biopsy_site <- factor(meta$biopsy_site,
                           levels = c("normal", "primary tumor", "colorectal cancer metastatic in the liver"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ biopsy_site)

# Run DESeq
dds <- DESeq(dds)

# Variance stabilizing transformation (for downstream analyses)
vsd <- vst(dds, blind = FALSE)



# Extract DESeq2 results (default: last level vs first level of biopsy_site)
res <- results(dds)

# Order by adjusted p-value (FDR)
res <- res[order(res$padj), ]

# Save DESeq2 results to CSV
write.csv(as.data.frame(res), file = "DESeq2_results.csv")

# Save normalized counts (VST-transformed)
norm_counts <- assay(vsd)
write.csv(norm_counts, file = "Normalized_Counts.csv")

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "biopsy_site", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(file.path(plot_dir, "PCA_plot.png"), width = 7, height = 6, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, color = biopsy_site)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples by Biopsy Site") +
  theme_minimal(base_size = 20) +   # bigger base font
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 11),
    legend.text  = element_text(size = 11),
    legend.title = element_blank(),      # remove legend title
    legend.position = "bottom"           # move legend to bottom
  )
dev.off()




# Save PCA data as CSV
write.csv(pcaData, file = "PCA_data.csv", row.names = TRUE)



# =====================
# Primary tumor vs Normal
# =====================
res_primary <- results(dds, contrast = c("biopsy_site", "primary tumor", "normal"))
res_primary <- lfcShrink(dds, coef = "biopsy_site_primary.tumor_vs_normal", type = "ashr", res = res_primary)
res_primary_ordered <- res_primary[order(res_primary$pvalue), ]
write.csv(as.data.frame(res_primary_ordered), "Primary_vs_Normal_DEGs.csv")

# MA Plot - Primary
png(file.path(plot_dir, "MA_plot_Primary_vs_Normal.png"),
    width = 15, height = 14, units = "in", res = 300)

# Set global text scaling
par(cex = 1.5)      # scales everything (axis text, labels, title, etc.)

# You can also fine-tune if needed:
# par(cex.axis = 1.5, cex.lab = 1.6, cex.main = 1.8)

plotMA(res_primary, ylim = c(-5, 5), 
       main = "MA Plot: Primary Tumor vs Normal")

dev.off()


# Volcano Plot - Primary
df_primary <- as.data.frame(res_primary)
df_primary <- df_primary[is.finite(df_primary$log2FoldChange) & is.finite(df_primary$padj), ]
df_primary$threshold <- ifelse(df_primary$padj < 0.05 & df_primary$log2FoldChange > 1, "Upregulated",
                               ifelse(df_primary$padj < 0.05 & df_primary$log2FoldChange < -1, "Downregulated", "NS"))

png(file.path(plot_dir, "Volcano_Primary_vs_Normal.png"), 
    width = 15, height = 14, units = "in", res = 300)

ggplot(df_primary, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7, size = 3) +   # bigger points
  scale_color_manual(values = c("Upregulated" = "green", 
                                "Downregulated" = "red", 
                                "NS" = "gray")) +
  geom_vline(xintercept = c(-1, 1), col = "blue", lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", lty = 2) +
  labs(title = "Volcano Plot: Primary Tumor vs Normal",
       x = "Log2Fold Change", 
       y = "-Log10 Adjusted P-value", 
       color = "DEG Status") +
  theme_minimal(base_size = 20) +   # consistent base font size
  theme(
    plot.title   = element_text(size = 24, hjust = 0.5), # big title
    axis.title   = element_text(size = 20),              # axis labels
    axis.text    = element_text(size = 16),              # axis ticks
    legend.title = element_text(size = 18),              # legend title
    legend.text  = element_text(size = 16),              # legend labels
    legend.position = "bottom"                           # put legend below
  )

dev.off()

  




# =====================
# Liver mets vs Normal
# =====================
res_liver <- results(dds, contrast = c("biopsy_site", 
                                       "colorectal cancer metastatic in the liver", 
                                       "normal"))
res_liver <- lfcShrink(dds, 
                       coef = "biopsy_site_colorectal.cancer.metastatic.in.the.liver_vs_normal", 
                       type = "ashr", res = res_liver)
res_liver_ordered <- res_liver[order(res_liver$pvalue), ]
write.csv(as.data.frame(res_liver_ordered), "LiverMets_vs_Normal_DEGs.csv")

# MA Plot - Liver
png(file.path(plot_dir, "MA_plot_Liver_vs_Normal.png"),
    width = 15, height = 14, units = "in", res = 300)

# Global text scaling for base R plot
par(cex = 1.5)

plotMA(res_liver, ylim = c(-5, 5), 
       main = "MA Plot: Liver Mets vs Normal")

dev.off()

# Volcano Plot - Liver
df_liver <- as.data.frame(res_liver)
df_liver <- df_liver[is.finite(df_liver$log2FoldChange) & is.finite(df_liver$padj), ]
df_liver$threshold <- ifelse(df_liver$padj < 0.05 & df_liver$log2FoldChange > 1, "Upregulated",
                             ifelse(df_liver$padj < 0.05 & df_liver$log2FoldChange < -1, "Downregulated", "NS"))

png(file.path(plot_dir, "Volcano_Liver_vs_Normal.png"), 
    width = 15, height = 14, units = "in", res = 300)

ggplot(df_liver, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7, size = 3) +   # bigger points
  scale_color_manual(values = c("Upregulated" = "green", 
                                "Downregulated" = "red", 
                                "NS" = "gray")) +
  geom_vline(xintercept = c(-1, 1), col = "blue", lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "blue", lty = 2) +
  labs(title = "Volcano Plot: Liver Mets vs Normal",
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       color = "DEG Status") +
  theme_minimal(base_size = 20) +
  theme(
    plot.title   = element_text(size = 24, hjust = 0.5),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.position = "bottom"
  )

dev.off()

# =====================
# Heatmap of top 20 variable genes
# =====================
# ============================
# Heatmaps: Top 20 variable genes
# - Unscaled (VST values)
# - Row Z-score scaled
# ============================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(grid)
})

# ---- 0) Parameters (tweak if needed) ----
n_top          <- 20
z_cap          <- 2.5     # cap z-scores at ±2.5 for robust visualization
dpi            <- 300
w_in           <- 8
h_in           <- 6
fs_main        <- 22
fs_label       <- 12
rowname_show   <- TRUE
clust_method   <- "ward.D2"   # good for heatmaps
dist_rows      <- "correlation"
dist_cols      <- "correlation"

# Diverging palette for Z-score; sequential for unscaled
pal_div  <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(101)
pal_cont <- colorRampPalette(brewer.pal(9, "YlGnBu"))(101)

# ---- 1) Build the data matrix (top variable genes) ----
# 'vsd' should be a SummarizedExperiment or DESeqTransform with assay(vsd)
mat_all <- assay(vsd)

# pick top variable genes (row variance)
sel_idx <- head(order(rowVars(mat_all), decreasing = TRUE), n_top)
heatmap_data <- as.matrix(mat_all[sel_idx, , drop = FALSE])

# drop any non-finite rows (rare but avoids plotting errors)
keep_rows <- apply(heatmap_data, 1, function(x) all(is.finite(x) & !is.na(x)))
heatmap_data <- heatmap_data[keep_rows, , drop = FALSE]

# ---- 2) Annotation (biopsy_site) ----
# Ensure annotation rows == columns of matrix
annotation_col <- data.frame(Biopsy = meta$biopsy_site)
rownames(annotation_col) <- colnames(heatmap_data)

# Nice, consistent colors for biopsy categories
biopsy_lvls <- unique(annotation_col$Biopsy)
biopsy_cols <- setNames(
  colorRampPalette(brewer.pal(max(3, min(9, length(biopsy_lvls))), "Set2"))(length(biopsy_lvls)),
  biopsy_lvls
)
ann_colors <- list(Biopsy = biopsy_cols)

# ---- 3) UN-SCALED (VST) HEATMAP ----
# Robust breaks using winsorized range (2nd–98th percentiles)
q <- quantile(heatmap_data, probs = c(0.02, 0.98), na.rm = TRUE)
vmin <- unname(q[1]); vmax <- unname(q[2])
if (vmin >= vmax) { # fallback if data are flat
  vmin <- min(heatmap_data, na.rm = TRUE)
  vmax <- max(heatmap_data, na.rm = TRUE)
}
breaks_cont <- seq(vmin, vmax, length.out = length(pal_cont) + 1)

png(file.path(plot_dir, "Heatmap_Top20Genes_unscaled.png"),
    width = 15, height = 14, units = "in", res = 300)
pheatmap(
  heatmap_data,
  annotation_col          = annotation_col,
  annotation_colors       = ann_colors,
  cluster_rows            = TRUE,
  cluster_cols            = TRUE,
  clustering_method       = clust_method,
  clustering_distance_rows= dist_rows,
  clustering_distance_cols= dist_cols,
  show_rownames           = rowname_show,
  fontsize                = fs_label,
  main                    = "Top Variable Genes (Unscaled VST)",
  color                   = pal_cont,
  breaks                  = breaks_cont,
  border_color            = NA,
  na_col                  = "grey90",
  angle_col               = 45
)
dev.off()

pheatmap(
  heatmap_data,
  annotation_col          = annotation_col,
  annotation_colors       = ann_colors,
  cluster_rows            = TRUE,
  cluster_cols            = TRUE,
  clustering_method       = clust_method,
  clustering_distance_rows= dist_rows,
  clustering_distance_cols= dist_cols,
  show_rownames           = rowname_show,
  fontsize                = fs_label,
  main                    = "Top Variable Genes (Unscaled VST)",
  color                   = pal_cont,
  breaks                  = breaks_cont,
  border_color            = NA,
  na_col                  = "grey90",
  angle_col               = 45
)
dev.off()

# ---- 4) Z-SCORE (row-scaled) HEATMAP ----
# Scale each gene (row) to mean 0, sd 1; cap extremes to ±z_cap for readability
z_mat <- t(scale(t(heatmap_data)))
# Replace any rows with sd=0 (NA after scaling) by zeros to retain structure
z_mat[!is.finite(z_mat)] <- 0
# cap values
z_mat[z_mat >  z_cap] <-  z_cap
z_mat[z_mat < -z_cap] <- -z_cap
breaks_div <- seq(-z_cap, z_cap, length.out = length(pal_div) + 1)

png(file.path(plot_dir, "Heatmap_Top20Genes_Zscore.png"),
    width = 15, height = 14, units = "in", res = 300)
pheatmap(
  z_mat,
  annotation_col          = annotation_col,
  annotation_colors       = ann_colors,
  cluster_rows            = TRUE,
  cluster_cols            = TRUE,
  clustering_method       = clust_method,
  clustering_distance_rows= dist_rows,
  clustering_distance_cols= dist_cols,
  show_rownames           = rowname_show,
  fontsize                = fs_label,
  main                    = "Top Variable Genes (Row Z-score)",
  color                   = pal_div,
  breaks                  = breaks_div,
  border_color            = NA,
  na_col                  = "grey90",
  angle_col               = 45
)
dev.off()


pheatmap(
  z_mat,
  annotation_col          = annotation_col,
  annotation_colors       = ann_colors,
  cluster_rows            = TRUE,
  cluster_cols            = TRUE,
  clustering_method       = clust_method,
  clustering_distance_rows= dist_rows,
  clustering_distance_cols= dist_cols,
  show_rownames           = rowname_show,
  fontsize                = fs_label,
  main                    = "Top Variable Genes (Row Z-score)",
  color                   = pal_div,
  breaks                  = breaks_div,
  border_color            = NA,
  na_col                  = "grey90",
  angle_col               = 45
)
dev.off()


# Save up/downregulated genes
up_primary <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange > 1)
down_primary <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_primary, "Upregulated_Primary_vs_Normal.csv")
write.csv(down_primary, "Downregulated_Primary_vs_Normal.csv")

up_liver <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange > 1)
down_liver <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_liver, "Upregulated_Liver_vs_Normal.csv")
write.csv(down_liver, "Downregulated_Liver_vs_Normal.csv")


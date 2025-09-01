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
setwd("C:/Users/newfaculty/Desktop/Bioinformatics_Club") ## Please change the working directory to where you have saved your meta and raw count csv data.
plot_dir <- "DESeq2_Plots"
dir.create(plot_dir, showWarnings = FALSE)

# Read input files
counts <- read.csv("raw_counts.csv", row.names = 1)  ## Check your raw count file name
meta   <- read.csv("MetaData.csv", row.names = 1)    ## Check your meta data file name

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
write.csv(counts(dds, normalized = TRUE), "Normalized_Counts.csv")

# Variance stabilizing transformation (for downstream analyses)
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
write.csv(norm_counts, file = "VSD_Counts.csv")


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
pcaData <- plotPCA(vsd, intgroup = "biopsy_site", ntop =1000, returnData = TRUE)
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

# Heatmaps: Top 20 Variable Genes
# Selects top 20 most variable genes (via row variance on VST matrix).
# Drop NAs, pick top 20 by smallest adjusted p-value
res_noNA <- res[!is.na(res$padj), ]
top_genes_tbl <- head(res_noNA[order(res_noNA$padj), , drop = FALSE], 20)
top_genes <- rownames(top_genes_tbl)

# Subset the expression matrix to those genes (using dds values)
hm_mat <- assay(vsd)[top_genes, , drop = FALSE]

# Heatmap
pheatmap(hm_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 20 DE genes (From Vst Normalized Gene Count)",
         filename = file.path(plot_dir, "VST_HeatmapTop20.png"),
         width = 8,
         height = 10,
         color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)
)

# Save up/downregulated genes
up_primary <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange > 1)
down_primary <- subset(res_primary_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_primary, "Upregulated_Primary_vs_Normal.csv")
write.csv(down_primary, "Downregulated_Primary_vs_Normal.csv")

up_liver <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange > 1)
down_liver <- subset(res_liver_ordered, padj < 0.05 & log2FoldChange < -1)
write.csv(up_liver, "Upregulated_Liver_vs_Normal.csv")
write.csv(down_liver, "Downregulated_Liver_vs_Normal.csv")









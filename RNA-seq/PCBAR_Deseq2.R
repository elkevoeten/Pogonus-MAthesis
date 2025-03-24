#### Differential gene expression analysis between P. chalceus ecotypes (Barbate) using Deseq2 ####

setwd("~/Documents/THESIS/RNA-seq-analysis")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ashr)

# Load count data
counts <- read.table("pogonus.counts.cleaned", header=TRUE, row.names=1, sep="\t")
counts <- counts[, -(1:5)]
head(counts)

# Define sample metadata
ecotype <- factor(c(rep("SW", 6), rep("LW", 6)))  # Adjust based on your samples
label <- data.frame(row.names=colnames(counts), condition=ecotype)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=label, design=~condition)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "LW", "SW"))  # LW vs SW

# Shrink log fold changes for visualization
resLFC <- lfcShrink(dds, coef="condition_SW_vs_LW", type="ashr")

# Save results
write.csv(as.data.frame(res), "DESeq2_results.csv")

# Sort results by log2 fold change
sorted_results <- res[order(res$log2FoldChange), ]

# Show top up- and down-regulated genes
head(sorted_results)  # Down-regulated genes
tail(sorted_results)  # Up-regulated genes

# PCA Plot
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  #geom_text(aes(label=rownames(pcaData)), vjust=-1, hjust=1, size=3) + # show labels > outlier is BAR4.01
  theme_minimal() +
  ggtitle("PCA of RNA-seq samples")

# MA Plot
plotMA(resLFC, main="MA plot", ylim=c(-5,5))

# Volcano Plot
res_df <- as.data.frame(res)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05)) +
  theme_minimal() +
  ggtitle("Volcano Plot")

# Heatmap of top differentially expressed genes
topGenes <- rownames(head(res[order(res$padj),], 30))
pheatmap(assay(rld)[topGenes, ], cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE)

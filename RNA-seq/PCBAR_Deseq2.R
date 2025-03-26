#### Differential gene expression analysis between P. chalceus ecotypes (Barbate) using Deseq2 ####

setwd("~/Documents/THESIS/RNA-seq-analysis")

library(DESeq2)
library(ggplot2)
library(ggnewscale)
library(pheatmap)
library(dplyr)
library(ashr)
library(RColorBrewer)
library(scales)
library(zoo)

#### DESEQ2 ####
# Load count data and only keep the count matrix
counts <- read.table("pogonus.counts.cleaned", header=TRUE, row.names=1, sep="\t")
positions <- counts[, (1:5)]
counts <- counts[, -(1:5)]
head(counts)

# Define sample metadata: SW is set as reference level
ecotype <- factor(c(rep("SW", 6), rep("LW", 6)))
ecotype
label <- data.frame(row.names=colnames(counts), condition=ecotype)
label

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=label, design=~condition)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "LW", "SW"))
res

# Amount of highly significantly DE genes
sum(res$padj < 0.01, na.rm=TRUE) #6881
sum(res$padj < 0.05, na.rm=TRUE) #8143
sum(res$padj < 0.001, na.rm=TRUE) #5588

# Save results
head(res)
write.table(res, file = "PCdeseq.results.tsv", sep = "\t", quote = FALSE, row.names = T)

# Merge positions and res table
positions$GeneID <- rownames(positions)
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
merged_res <- merge(res_df, positions, by="GeneID", all.x=TRUE)
sum(merged_res$padj < 0.01, na.rm=TRUE) #6881

# Filter for padj < 0.01
sig_res <- merged_res[merged_res$padj < 0.01, ]
sum(sig_res$padj < 0.01, na.rm=TRUE) #6881
length(unique(sig_res$GeneID)) #6882
nrow(sig_res) #7076

# Check why there's a difference in rows
table(duplicated(positions$GeneID))  # Check for duplicates in positions
table(duplicated(res_df$GeneID))  # Check for duplicates in res_df
table(duplicated(merged_res$GeneID))  # Check for duplicates after merging
sum(is.na(merged_res$padj)) #195
sum(is.na(sig_res$padj)) #195 => it's the NAs! 

sig_res <- merged_res[!is.na(merged_res$padj) & merged_res$padj < 0.01, ]
nrow(sig_res) #6881
length(unique(sig_res$GeneID)) #6881

head(sig_res)

# Filter out the genes that have -1 < LFC < 1
LFC_sig_res <- sig_res[sig_res$log2FoldChange < -1 | sig_res$log2FoldChange > 1, ]
nrow(LFC_sig_res) #4324
write.table(LFC_sig_res, file = "PCdeseq.LFC_sig_res.tsv", sep = "\t", quote = FALSE, row.names = T)


# Make the plot to show the LFC against the position in the genome
summary(-log10(LFC_sig_res$padj))

ggplot(LFC_sig_res, aes(x = Start, y = log2FoldChange)) +
  geom_point(aes(color = log2FoldChange > 0, size = -log10(padj)), alpha = 0.8) +
  scale_color_manual(
    values = c("FALSE" = "cornflowerblue", "TRUE" = "brown3"),
    labels = c("FALSE" = "Downregulated", "TRUE" = "Upregulated")
    ) +
  scale_size_continuous(
    range = c(0.1, 6), 
    breaks = c(1000, 100, 50, 25, 10), 
    labels = c("1000", "100", "50", "25", "10"),
    name = "-log10(padj)"
    ) +
  facet_wrap(~Chr, scales = "free") + 
  theme_minimal() +
  labs(
    title = "Log2 Fold Change vs. Gene Start Position",
    x = "Gene position (1e7 bp)",
    y = "Log2 Fold Change",
    color = "Log2FoldChange",
    size = "-log10(padj)"
    ) +
  scale_x_continuous(labels = function(x) x / 1e7) +
  theme(legend.position = "bottom")

ggsave("LFC_sig_plot.pdf", width = 15, height = 9)


#### PLOT FST SIGNALS AS BACKGROUND ####
fst_dir <- "/Users/elkevoeten/Documents/THESIS/RNA-seq-analysis/fst-signals"
fst_files <- list.files(fst_dir, pattern = "*.stats", full.names = TRUE)
fst_data_list <- lapply(fst_files, read.csv)
fst_data_combined <- do.call(rbind, fst_data_list)
fst_data_combined <- fst_data_combined %>%
  rename(Chr = scaffold, Start = start, Fst = Fst_S_HUE_T_S_COT_S) %>%
  filter(Fst >= 0)  # Keep only non-negative Fst values


# Fst plot trial
ggplot(fst_data_combined, aes(x = Start, y = Fst, color=Fst)) +
  geom_point(alpha=0.8) +
  scale_color_gradient(low = "grey", high = "grey24") + 
  facet_wrap(~Chr, scales = "free_x") + 
  labs(
    title = "Fst Across Chromosomes",
    x = "Gene position",
    y = "Fst"
    ) +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Rotate and adjust side labels
    legend.position = "bottom"
    ) +
  scale_x_continuous(labels = function(x) x / 1e7)


# Combine the plots
ggplot() +
  geom_point(data=fst_data_combined, aes(x= Start, y=Fst), color = "grey", alpha = 0.5) +
  geom_point(data=LFC_sig_res, aes(x=Start, y=log2FoldChange/20, color = log2FoldChange > 0, size = -log10(padj)), alpha = 0.8) +
  scale_color_manual(
    values = c("FALSE" = "cornflowerblue", "TRUE" = "brown3"),
    labels = c("FALSE" = "Downregulated", "TRUE" = "Upregulated")
    ) +
  scale_size_continuous(
    range = c(0.1, 6), 
    breaks = c(1000, 100, 50, 25, 10), 
    labels = c("1000", "100", "50", "25", "10"),
    name = "-log10(padj)"
    ) +
  scale_y_continuous(
    sec.axis = sec_axis(trans = ~.*20, name = "log2FoldChange")
    ) +
  scale_x_continuous(labels = function(x) x / 1e7) +
  facet_wrap(~Chr, scales = "free_x") +
  labs(
    title = "Fst and Log2FoldChange plot",
    x = "Position (10^7 bp)",
    color = "Log2FoldChange",
    size = "-log10(padj)"
    ) +
  theme_minimal()

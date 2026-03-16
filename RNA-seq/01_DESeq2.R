cat("\014")   #clear console
rm(list=ls()) #clear environment

library(rstudioapi)
main_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(main_path))

intermediate_output_path <- "intermediate-output"
if (!dir.exists(intermediate_output_path)) {
  dir.create(intermediate_output_path)
  cat("Output path created:", intermediate_output_path, "\n")
} else {
  cat("Output path already exists:", intermediate_output_path, "\n")
}

############################################################################################################
# script used to run differential gene expression analysis between Pogonus ecotypes (Barbate) using DESeq2 #
# note: the counts matrix is the output of the upstream analysis that aligns the reads to the reference    #
# genome and uses featureCounts to count the amount of expressed genes. DESeq2 quantifies the differential #
# expressed genes. I rewrote the original script so only all samples (ALL) and the P. Chalceus and P.      #
# Gilvipes (CG) comparisons are still included.                                                            #
############################################################################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")

library(DESeq2)       #for DESeq2 analysis
library(apeglm)       #for lfcShrink calculations
library(ggplot2)
library(dplyr)        #for dataframe manipulation
library(pheatmap)     #for heatmap plotting

# DATA LOADING ---------------------------------------------------------------------------------------------
counts <- read.table("../data/pogonus.counts.cleaned", header=TRUE, row.names=1, sep="\t")
positions <- counts[, (1:5)] # extract gene position info
positions$GeneID <- rownames(positions) # get geneIDs for later merging
counts <- counts[, -(1:5)]   # get count matrix for downstream analysis

# DESEQ2 ---------------------------------------------------------------------------------------------------
# [A] ADD METADATA -----------------------------------------------------------------------------------------
ALL.counts <- counts
ALL.metadata <- data.frame(
  row.names = colnames(ALL.counts),
  ecotype = factor(c(rep("SW", 6), rep("LW", 6))),
  habitat = factor(c(rep("tidal", 6), rep("seasonal", 6))),
  species = factor(c(rep("chalceus", 6), "littoralis", rep("gilvipes", 5)))
  )

CG.counts <- counts[, c(1:6, 8:12)]
CG.metadata <- data.frame(
  row.names = colnames(CG.counts),
  species = factor(c(rep("chalceus", 6), rep("gilvipes", 5)))
  )

# [B] CREATE DESEQ DATASET ---------------------------------------------------------------------------------
# ALL SAMPLES
ALL.dds <- DESeqDataSetFromMatrix(countData=ALL.counts, colData=ALL.metadata, design=~habitat)
ALL.dds$habitat <- factor(ALL.dds$habitat, levels = c("seasonal","tidal")) #set reference level to 'seasonal'
ALL.dds <- ALL.dds[rowSums(counts(ALL.dds) > 10) >= 6, ] #pre-filter low-count genes

# CHALCEUS VS. GILVIPES
CG.dds <- DESeqDataSetFromMatrix(countData=CG.counts, colData=CG.metadata, design=~species)
CG.dds$species <- factor(CG.dds$species, levels = c("gilvipes", "chalceus")) #set reference level to 'seasonal species'
CG.dds <- CG.dds[rowSums(counts(CG.dds) > 10) >= 5, ] #pre-filter low-count genes

# [C] RUN DESEQ2  ------------------------------------------------------------------------------------------
# ALL SAMPLES
ALL.dds <- DESeq(ALL.dds)
ALL.res <- results(ALL.dds)
ALL.res

ggplot(as.data.frame(ALL.res), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = padj < 0.01), size = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(name = "padj < 0.01", values = setNames(c("red", "black"), c(T, F)))

ggplot(as.data.frame(ALL.res), aes(x = baseMean, y = lfcSE)) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = "log10") +
  ylim(0, 3)

resultsNames(ALL.dds)
ALL.resLFC <- lfcShrink(ALL.dds, coef="habitat_tidal_vs_seasonal")
ALL.resLFC

ggplot(as.data.frame(ALL.resLFC), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = padj < 0.01), size = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(name = "padj < 0.01", values = setNames(c("red", "black"), c(T, F)))

ggplot(as.data.frame(ALL.resLFC), aes(x = baseMean, y = lfcSE)) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = "log10") +
  ylim(0, 3)

# CHALCEUS VS. GILVIPES
CG.dds <- DESeq(CG.dds)
CG.res <- results(CG.dds)
CG.res

ggplot(as.data.frame(CG.res), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = padj < 0.01), size = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(name = "padj < 0.01", values = setNames(c("red", "black"), c(T, F)))

ggplot(as.data.frame(CG.res), aes(x = baseMean, y = lfcSE)) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = "log10") +
  ylim(0, 3)

resultsNames(CG.dds)
CG.resLFC <- lfcShrink(CG.dds, coef="species_chalceus_vs_gilvipes")
CG.resLFC

ggplot(as.data.frame(CG.resLFC), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = padj < 0.01), size = 0.5) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(name = "padj < 0.01", values = setNames(c("red", "black"), c(T, F)))

ggplot(as.data.frame(CG.resLFC), aes(x = baseMean, y = lfcSE)) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = "log10") +
  ylim(0, 3)

# DOWNSTREAM ANALYSIS ---------------------------------------------------------------------------------------
# ALL SAMPLES
# save results table
head(ALL.res)
# write.table(ALL.res, file = "intermediate-output/DESeq2.ALL.results.tsv", sep = "\t", quote = FALSE, row.names = T)
head(ALL.resLFC)
# write.table(ALL.resLFC, file = "intermediate-output/DESeq2.ALL.resultsLFC.tsv", sep = "\t", quote = FALSE, row.names = T)

ALL.resLFC.df <- as.data.frame(ALL.resLFC)
ALL.resLFC.df$GeneID <- rownames(ALL.resLFC.df)
ALL.resLFC.df <- merge(ALL.resLFC.df, positions, by="GeneID", all.x=TRUE)
sum(ALL.resLFC.df$padj < 0.01, na.rm=TRUE) #6398
ALL.resLFC.sig <- ALL.resLFC.df[!is.na(ALL.resLFC.df$padj) & ALL.resLFC.df$padj < 0.01, ] # Filter for padj < 0.01
nrow(ALL.resLFC.sig) #6398

ALL.resLFC.sig.DEG <- ALL.resLFC.sig[ALL.resLFC.sig$log2FoldChange < -1 | ALL.resLFC.sig$log2FoldChange > 1, ]
nrow(ALL.resLFC.sig.DEG) #3719
sum(ALL.resLFC.sig.DEG$log2FoldChange > 1) #2793
sum(ALL.resLFC.sig.DEG$log2FoldChange < -1) #926
# write.table(ALL.resLFC.sig.DEG, file = "intermediate-output/DESeq2.ALL.resLFC.sig.DEG.tsv", sep = "\t", quote = FALSE, row.names = T)

# plot heatmap
ALL.upregulated_genes <- ALL.resLFC.df %>%
  filter(log2FoldChange > 1, padj < 0.01) %>%
  pull(GeneID)
ALL.downregulated_genes <- ALL.resLFC.df %>% 
  filter(log2FoldChange < 1, padj < 0.01) %>%
  pull(GeneID)
ALL.resLFC.df <- ALL.resLFC.df %>% arrange(log2FoldChange) 

pheatmap(counts[c(ALL.upregulated_genes, ALL.downregulated_genes), ], scale = "row", show_rownames = F)

# CHALCEUS VS. GILVIPES
# save results table
head(CG.res)
# write.table(CG.res, file = "intermediate-output/DESeq2.CG.results.tsv", sep = "\t", quote = FALSE, row.names = T)
head(CG.resLFC)
# write.table(CG.resLFC, file = "intermediate-output/DESeq2.CG.resultsLFC.tsv", sep = "\t", quote = FALSE, row.names = T)

CG.resLFC.df <- as.data.frame(CG.resLFC)
CG.resLFC.df$GeneID <- rownames(CG.resLFC.df)
CG.resLFC.df <- merge(CG.resLFC.df, positions, by="GeneID", all.x=TRUE)
sum(CG.resLFC.df$padj < 0.01, na.rm=TRUE) #7380 
CG.resLFC.sig <- CG.resLFC.df[!is.na(CG.resLFC.df$padj) & CG.resLFC.df$padj < 0.01, ] # Filter for padj < 0.01
nrow(CG.resLFC.sig) #7380

CG.resLFC.sig.DEG <- CG.resLFC.sig[CG.resLFC.sig$log2FoldChange < -1 | CG.resLFC.sig$log2FoldChange > 1, ]
nrow(CG.resLFC.sig.DEG) #4315
sum(CG.resLFC.sig.DEG$log2FoldChange > 1) #3211
sum(CG.resLFC.sig.DEG$log2FoldChange < -1) #1104
# write.table(CG.resLFC.sig.DEG, file = "intermediate-output/DESeq2.CG.resLFC.sig.DEG.tsv", sep = "\t", quote = FALSE, row.names = T)

# plot heatmap
CG.upregulated_genes <- CG.resLFC.df %>%
  filter(log2FoldChange > 1, padj < 0.01) %>%
  pull(GeneID)
CG.downregulated_genes <- CG.resLFC.df %>% 
  filter(log2FoldChange < 1, padj < 0.01) %>%
  pull(GeneID)
CG.resLFC.df <- CG.resLFC.df %>% arrange(log2FoldChange) 

pheatmap(counts[c(CG.upregulated_genes, CG.downregulated_genes), ], scale = "row", show_rownames = F)

# MA Plot
plotMA(CG.resLFC, main="MA plot", ylim=c(-5,5))

# Volcano Plot
ggplot(CG.resLFC, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05)) +
  theme_minimal() +
  ggtitle("Volcano Plot")

library(ggplot2)

output_path <- "output-figures"
if (!dir.exists(output_path)) {
  dir.create(output_path)
  cat("Output path created:", output_path, "\n")
} else {
  cat("Output path already exists:", output_path, "\n")
}

############################################################################################################
# script to run a PCA analysis for all samples of the RNA seq analysis; uses ALL.dds object from 01_DESeq2.R
############################################################################################################

# PCA CALCULATIONS -----------------------------------------------------------------------------------------
# Extract rlog-transformed expression matrix
ALL.rld <- rlog(ALL.dds, blind=TRUE)
expr_matrix <- assay(ALL.rld)  # genes (rows) × samples (columns)
expr_matrix

# Transpose to samples × genes for PCA
pca <- prcomp(t(expr_matrix), center = TRUE, scale. = FALSE)
pca

summary(pca) # Calculate % variance explained
percentVar <- round(100 * summary(pca)$importance[2, 1:2], 1)
percentVar
plot(pca, type = "l", main = "Scree Plot of PCA")

# PCA scores (per sample; used for scatter plots)
pca_scores <- as.data.frame(pca$x)  # rows = samples, columns = PCs
head(pca_scores)

# PCA loadings (eigenvectors; gene contributions to each PC)
pca_loadings <- as.data.frame(pca$rotation)  # rows = genes, columns = PCs
head(pca_loadings)

# PCA-PLOT ----------------------------------------------------------------------------------------------------
pca_scores$sample <- colnames(ALL.rld)
pca_scores$species <- c(rep("P. chalceus", 6), "P. littoralis", rep("P. gilvipes", 5))
pca_scores$habitat <- c(rep("tidal", 6), rep("seasonal", 6))
pca_scores

PCA.plot <- ggplot(pca_scores, aes(PC1, PC2, color = habitat)) +
  geom_point(size = 3, alpha = 0.95) +
  # geom_text(aes(label = rownames(pca_scores)), vjust = -1, hjust = 1, size = 3) + 
  theme_minimal() +
  ggtitle("PCA of RNA-seq samples") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  ) +
  scale_color_manual(
    values = c("seasonal" = "#E97132", "tidal" = "#4E95D9")
  ) +
  guides(color = guide_legend(title = "Habitat"))
PCA.plot

# ggsave(PCA.plot,
#        filename = "output-figures/PCA.plot.pdf",
#        device = "pdf", 
#        width = 20, height = 14, units = "cm")

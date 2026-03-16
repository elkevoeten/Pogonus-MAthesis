cat("\014")   #clear console
rm(list=ls()) #clear environment

library(rstudioapi)
main_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(main_path))

output_path <- "output-figures"
if (!dir.exists(output_path)) {
  dir.create(output_path)
  cat("Output path created:", output_path, "\n")
} else {
  cat("Output path already exists:", output_path, "\n")
}

############################################################################################################
# script to visualise the differential expressed genes (Chalceus vs. Gilvipes)
############################################################################################################

CG.resLFC.sig.DEG.plot <- read.delim("intermediate-output/DESeq2.CG.resLFC.sig.DEG.annotated.tsv", sep="\t")
inversions <- read.delim("../data/crossings_inversions.tsv", header = TRUE, sep = "\t")
regions <- read.delim("../data/pogonus_regions.tsv", header = TRUE, sep = "\t")
fst_intervals <- read.delim("intermediate-output/fst-intervals.tsv", sep = "\t")
fst_data_combined <- read.delim("intermediate-output/fst-data.tsv", sep = "\t")
repeat_content <- read.delim("intermediate-output/repeat-content.tsv", sep = "\t")

library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(glue)
library(cowplot)

# PLOT DEGS TO GENOME POSITION (WITH REPEAT CONTENT AND FST SIGNALS) ---------------------------------------
# List of chromosomes to loop through
chromosomes <- unique(CG.resLFC.sig.DEG.plot$chromosome)

# Store plots and widths
plot_list <- list()
widths_vec <- numeric(length(chromosomes))

# Function to conditionally set y-axis text and ticks
get_yaxis_label <- function(chr) {
  if (chr == "Chr1") element_text() else element_blank()
}
get_yaxis_text <- function(chr) {
  if (chr == "Chr1") element_text() else element_blank()
}
get_yaxis_ticks <- function(chr) {
  if (chr == "Chr1") element_line() else element_blank()
}

for (i in seq_along(chromosomes)) {
  chr <- chromosomes[i]
  chr_label <- chr
  chr_safe <- gsub("\\.", "_", chr_label)
  
  chr_inversions <- inversions %>% filter(chromosome == chr)
  chr_data_up <- CG.resLFC.sig.DEG.plot %>% filter(chromosome == chr, diff_expr == "upregulation")
  chr_data_down <- CG.resLFC.sig.DEG.plot %>% filter(chromosome == chr, diff_expr == "downregulation")
  chr_fst <- fst_data_combined %>% filter(chromosome == chr)
  chr_repeat_content <- repeat_content %>% filter(chromosome == chr)
  
  xlims_chr <- range(c(chr_data_up$Start, chr_data_down$Start), na.rm = TRUE)
  chr_width <- diff(xlims_chr)
  widths_vec[i] <- chr_width
  
  # Plot 1: Upregulated genes
  p1_up <- ggplot(chr_data_up) +
    geom_rect(data = chr_inversions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#ADE0FF", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(aes(x = Start, y = log2FoldChange, color = location, size = -log10(padj)), alpha = 0.5, show.legend = TRUE) +
    geom_boxplot(
      data = chr_data_up %>% group_by(region) %>% mutate(region_pos = median(Start, na.rm = TRUE)),
      aes(x = region_pos, y = log2FoldChange, group = region, color = location),
      alpha = 0.9, position = position_nudge(x = 0), width = 5e6, outlier.shape = NA, show.legend = FALSE
    ) +
    labs(title = glue("{chr_label}")) +
    scale_x_continuous(labels = function(x) x / 1e7) +
    coord_cartesian(xlim = xlims_chr) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = get_yaxis_label(chr),
      axis.text.y = get_yaxis_text(chr),
      axis.ticks.y = get_yaxis_ticks(chr),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#F5FBFF", color = NA)
    ) +
    ylim(0, 16)
  
  
  # Plot 2: Downregulated genes
  p1_down <- ggplot(chr_data_down) +
    geom_rect(data = chr_inversions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#ADE0FF", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(aes(x = Start, y = log2FoldChange, color = location, size = -log10(padj)), alpha = 0.5) +
    geom_boxplot(
      data = chr_data_down %>% group_by(region) %>% mutate(region_pos = median(Start, na.rm = TRUE)),
      aes(x = region_pos, y = log2FoldChange, group = region, color = location), 
      alpha = 0.9, position = position_nudge(x = 0), width = 5e6, outlier.shape = NA, show.legend = FALSE
    ) +
    scale_x_continuous(labels = function(x) x / 1e7) +
    coord_cartesian(xlim = xlims_chr) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = get_yaxis_label(chr),
      axis.text.y = get_yaxis_text(chr),
      axis.ticks.y = get_yaxis_ticks(chr),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#F5FBFF", color = NA)
    ) +
    ylim(-8, 0)  # Set y-axis range for downregulated genes
  
  # Plot 3: Repeat content
  p2 <- ggplot(chr_repeat_content, aes(x = start, y = perc)) +
    geom_rect(data = chr_inversions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#ADE0FF", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(color = "orange4", size = 1.5, alpha = 0.7) +
    labs(x = "", y = "Repeat %") +
    scale_x_continuous(labels = function(x) x / 1e7) +
    coord_cartesian(xlim = xlims_chr) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = get_yaxis_label(chr),
      axis.text.y = get_yaxis_text(chr),
      axis.ticks.y = get_yaxis_ticks(chr),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#F5FBFF", color = NA)
    )
  
  # Plot 4: FST
  p3 <- ggplot(chr_fst, aes(x = Start, y = Fst, color = Fst)) +
    geom_rect(data = chr_inversions, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#ADE0FF", alpha = 0.3, inherit.aes = FALSE) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = 0.2, linetype = "dashed") +
    scale_color_gradient(low = "#ADDDFF", high = "#00477A", guide = "none") +
    scale_x_continuous(labels = function(x) x / 1e7) +
    coord_cartesian(xlim = xlims_chr) +
    labs(x = expression("Position ("*10^7*" bp)"), y = "FST") +
    theme_minimal() +
    theme(
      axis.title.y = get_yaxis_label(chr),
      axis.text.y = get_yaxis_text(chr),
      axis.ticks.y = get_yaxis_ticks(chr),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      panel.background = element_rect(fill = "#F5FBFF", color = NA)
    ) +
    ylim(0, 1)  # Set y-axis range for FST plot
  
  # Combine vertically per chromosome
  chr_plot <- (p1_up / p1_down / p2 / p3) + plot_layout(ncol = 1, heights = c(3.5, 3.5, 1.5, 1.5))
  
  plot_list[[i]] <- chr_plot
}

# Normalize widths
widths_vec <- widths_vec / max(widths_vec)

# Combine all chromosomes horizontally with proportional widths
final_plot <- wrap_plots(plot_list, nrow = 1, widths = widths_vec)

# make combined x-label
x_label <- wrap_elements(
  grid::textGrob(
    expression("Position ("*10^7*" bp)"),
    gp = grid::gpar(fontsize = 14),
    just = "center"
  )
)

sample_plot <- ggplot(data = chr_data_up) +
  geom_point(aes(x = Start, y = log2FoldChange, color = location, size = -log10(padj)), alpha = 0.5) +
  labs(title = "Legend for Combined Plots") +
  theme_minimal() +
  theme(legend.position = "top")

# Extract the legend
legend_plot <- get_legend(sample_plot + theme(legend.text = element_text(size = 14),
                                              legend.title = element_text(size = 16)))

# Create a blank ggplot for the legend part
legend_ggplot <- ggplot() + 
  theme_void() +  # This clears any default settings
  annotation_custom(as_grob(legend_plot), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# Combine the blank legend with the final plot and x_label
final_combined_plot <- (legend_ggplot / final_plot / x_label) +
  plot_layout(heights = c(1, 10, 0.5))  # Adjust heights if needed
final_combined_plot

# Save the final plot
# ggsave("output-figures/DEGs.plot.pdf",
#        plot = final_combined_plot,
#        device = "pdf",
#        width = 60, height = 25, units = "cm")

# WHOLE GENOME BOXPLOT -------------------------------------------------------------------------------------
# Create individual plots
p_up <- ggplot(data = CG.resLFC.sig.DEG.plot %>% filter(diff_expr == "upregulation"),
               aes(x = location, y = log2FoldChange, group = location, colour = location)) +
  geom_boxplot(alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 17, label.x = 1.4) +
  coord_cartesian(ylim = c(0, 18)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = c(
    "background" = "Background",
    "inversion" = "Putative inversion")) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_down <- ggplot(data = CG.resLFC.sig.DEG.plot %>% filter(diff_expr == "downregulation"),
                 aes(x = location, y = log2FoldChange, group = location, colour = location)) +
  geom_boxplot(alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = -10, label.x = 1.4) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = c(
    "background" = "Background",
    "inversion" = "Putative inversion")) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Stack the plots
CG.boxplots <- p_up / p_down
CG.boxplots

# ggsave("output-figures/DEGs.boxplot.pdf",
#        plot = CG.boxplots,
#        device = "pdf",
#        width = 4, height = 8, units = "cm")

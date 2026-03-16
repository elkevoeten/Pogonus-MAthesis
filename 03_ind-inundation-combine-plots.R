library(patchwork)

# Combine the two plots side-by-side
SSinundation_combined <- ind_inundation_september + ind_inundation_apr_species +
  plot_layout(ncol = 2, widths = c(3, 5))  # 2 columns: September on left, April on right

# Display the combined plot
SSinundation_combined

# Save the combined figure
ggsave("ind_inundation_combined.jpeg", plot = SSinundation_combined, width = 32, height = 12, units = "cm", dpi = 600)



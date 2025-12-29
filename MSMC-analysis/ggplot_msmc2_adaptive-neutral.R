library(here)
base_dir <- here("/Users/elkevoeten/Library/CloudStorage/OneDrive-Personal/Documenten/Biology-KUL/THESIS/msmc/ggplot-msmc-adaptive-neutral")

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# define populations
pops <- list(
  Belgium = c(
    'GC129388', 'GC129389', 'GC129390', 'GC129391', 'GC129392', 'GC129393', 'GC129394', 'GC129395', 'GC129396', 'GC129397', 'GC129398', 'GC129399',
    'GC136084', 'GC136085', 'GC136086', 'GC136087', 'GC136088', 'GC136089', 'GC136090', 'GC136091', 'GC136092', 'GC136093', 'GC136094', 'GC136095'),
  France = c(
    'GC129400', 'GC129401', 'GC129402', 'GC129403', 'GC129404', 'GC129405', 'GC129406', 'GC129407', 'GC129408', 'GC129409', 'GC129410', 'GC129411',
    'GC129422', 'GC129423', 'GC136078', 'GC136079', 'GC136080', 'GC136081', 'GC136082', 'GC136083', 'GC136096', 'GC136097', 'GC136098', 'GC136099', 
    'GC136100', 'GC136101'),
  Portugal = c(
    'GC129412', 'GC129413', 'GC129414', 'GC129415', 'GC129416', 'GC129417', 'GC129418', 'GC129419', 'GC129420', 'GC129421'),
  Spain = c(
    'GC136107', 'GC136108', 'GC136109', 'GC136110', 'GC136111', 'GC136112', 'GC136113', 'GC136114', 'GC136115', 'GC136116'),
  Heist = c(
    'GC129427', 'GC129428', 'GC129429', 'GC129430', 'GC129431', 'GC129432', 'GC129433', 'GC129434', 'GC129435', 'GC129437', 'GC129438', 'GC136117', 
    'GC136118', 'GC136119', 'GC136120', 'GC136121', 'GC136122', 'GC136123', 'GC136124', 'GC136125', 'GC136126', 'GC136127', 'GC136128'),
  UK = c(
    'GC129424', 'GC129425', 'GC129426', 'GC136105', 'GC136106'),
  Mediterranean = c(
    'GC136102', 'GC136103', 'GC136104')
)

tidal = c(
  'GC129388', 'GC129389', 'GC129390', 'GC129391', 'GC129392', 'GC129393', 'GC129400', 'GC129401', 'GC129402', 'GC129403', 'GC129404', 'GC129405', 
  'GC129412', 'GC129413', 'GC129414', 'GC129415', 'GC129416', 'GC129424', 'GC129425', 'GC129426', 'GC136090', 'GC136091', 'GC136092', 'GC136093', 
  'GC136094', 'GC136095', 'GC136105', 'GC136106', 'GC136107', 'GC136108', 'GC136114', 'GC136115', 'GC136116', 'GC136078', 'GC136079', 'GC136080', 
  'GC136081', 'GC136082', 'GC136083', 'GC136123', 'GC136124', 'GC136125', 'GC136126', 'GC136127', 'GC136128', 'GC129433', 'GC129434', 'GC129435', 
  'GC129437', 'GC129438')
seasonal = c(
  'GC129394', 'GC129395', 'GC129396', 'GC129397', 'GC129398', 'GC129399', 'GC129406', 'GC129407', 'GC129408', 'GC129409', 'GC129410', 'GC129411', 
  'GC129417', 'GC129418', 'GC129419', 'GC129420', 'GC129421', 'GC129422', 'GC129423', 'GC136084', 'GC136085', 'GC136086', 'GC136087', 'GC136088', 
  'GC136089', 'GC136096', 'GC136097', 'GC136098', 'GC136099', 'GC136100', 'GC136101', 'GC136102', 'GC136103', 'GC136104', 'GC136109', 'GC136110', 
  'GC136111', 'GC136112', 'GC136113', 'GC129427', 'GC129428', 'GC129429', 'GC129430', 'GC129431', 'GC129432', 'GC136117', 'GC136118', 'GC136119', 
  'GC136120', 'GC136121','GC136122')

sample_to_ecotype <- function(sample) {
  case_when(
    sample %in% tidal    ~ "Tidal",
    sample %in% seasonal ~ "Seasonal",
    TRUE ~ NA_character_
  )
}

# LOAD DATA
adaptive_dir <- file.path(base_dir, "msmc-adaptive")
neutral_dir  <- file.path(base_dir, "msmc-neutral")

read_msmc <- function(file, sample, pop, region) {
  df <- read.table(file, header = TRUE)
  df %>%
    transmute(
      time     = leftX,
      Ne       = popSize,
      sample   = sample,
      popGroup = pop,
      ecotype  = sample_to_ecotype(sample),
      region   = region
    )
}

# loop over populations, samples and regions
plot_data_list <- list()

for (pop in names(pops)) {
  for (sample in pops[[pop]]) {
    # adaptive regions
    adaptive_file <- file.path(adaptive_dir, paste0(sample, ".final.Rin.txt"))
    if (file.exists(adaptive_file)) {
      plot_data_list[[length(plot_data_list) + 1]] <-
        read_msmc(adaptive_file, sample, pop, "Adaptive")
    }
    # neutral regions
    neutral_file <- file.path(neutral_dir, paste0(sample, ".final.Rin.txt"))
    if (file.exists(neutral_file)) {
      plot_data_list[[length(plot_data_list) + 1]] <-
        read_msmc(neutral_file, sample, pop, "Neutral")
    }
  }
}

# combine everything into one dataframe
plot_data <- bind_rows(plot_data_list)

plot_data <- plot_data %>%
  filter(
    time > 0,
    Ne > 0,
    !is.na(ecotype)
  )

plot_data <- plot_data %>%
  filter(!popGroup %in% c("UK", "Mediterranean"))

# --------------------------------------------------------------------------------------------------------------------
library(cowplot)

msmc_plot_adaptive <- function(df) {
  ggplot(df, aes(x = time, y = Ne)) +
    geom_step(aes(group = sample, color = ecotype), alpha = 0.3, linewidth = 0.5) +
    geom_smooth(aes(group = ecotype, color = ecotype), method = "loess", se = FALSE, linewidth = 1.4, span = 0.5) +
    geom_vline(xintercept = 190000, colour = "darkgrey", linewidth = 0.5) +
    scale_x_log10(
      breaks = 10^(2:7),
      labels = function(x) parse(text = paste0("10^", log10(x)))
    ) +
    scale_y_log10(
      breaks = 10^(2:8),
      labels = function(x) parse(text = paste0("10^", log10(x)))
    ) +
    scale_color_manual(
      values = c("Tidal" = "#4E95D9", "Seasonal" = "#E97132"),
      labels = c("LW", "SW")
    ) +
    facet_wrap(~ popGroup, scales = "free_y", nrow = 1) +
    labs(x = NULL, y = expression(N[e])) +  # geen x-as titel
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}

msmc_plot_neutral <- function(df) {
  ggplot(df, aes(x = time, y = Ne)) +
    geom_step(aes(group = sample, color = ecotype), alpha = 0.3, linewidth = 0.5) +
    geom_smooth(aes(group = ecotype, color = ecotype), method = "loess", se = FALSE, linewidth = 1.4, span = 0.5) +
    geom_vline(xintercept = 190000, colour = "darkgrey", linewidth = 0.5) +
    scale_x_log10(
      breaks = 10^(2:7),
      labels = function(x) parse(text = paste0("10^", log10(x)))
    ) +
    scale_y_log10(
      breaks = 10^(2:8),
      labels = function(x) parse(text = paste0("10^", log10(x)))
    ) +
    scale_color_manual(
      values = c("Tidal" = "#4E95D9", "Seasonal" = "#E97132"),
      labels = c("LW", "SW")
    ) +
    facet_wrap(~ popGroup, scales = "free_y", nrow = 1) +
    labs(
      x = expression(paste("Years ago (g = 1, ", mu, " = 2.1 × 10"^{-9}, ")")),
      y = expression(N[e]),
      color = "Ecotype"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 18),
      legend.title = element_text(face = "bold", size = 18),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid = element_blank()
    )
}

library(grid)

right_title_adaptive <- textGrob("Adaptive regions", rot = -90, gp = gpar(fontsize = 18, fontface = "bold"))
right_title_neutral <- textGrob("Neutral regions", rot = -90, gp = gpar(fontsize = 18, fontface = "bold"))

adaptive_plot <- msmc_plot_adaptive(
  plot_data %>% filter(region == "Adaptive")
) + inset_element(right_title_adaptive, left = 1, right = 1.07, bottom = 0, top = 1)

neutral_plot <- msmc_plot_neutral(
  plot_data %>% filter(region == "Neutral")
) + inset_element(right_title_neutral, left = 1, right = 1.07, bottom = 0, top = 1)

library(patchwork)

combined_plot_AN <- adaptive_plot / 
  plot_spacer() / 
  neutral_plot +
  plot_layout(heights = c(1, 0.05, 1))
combined_plot_AN

ggsave("msmc_adaptive-neutral.pdf", plot = combined_plot_AN, width = 12, height = 7)
ggsave("msmc_adaptive-neutral.jpeg", plot = combined_plot_AN, width = 12, height = 7, dpi = 600)

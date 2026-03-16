cat("\014")   #clear console
rm(list=ls()) #clear environment

library(rstudioapi)
main_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(main_path))
getwd()

############################################################################################################
# script to extract points in the genome where Fst is elevated to later differentiate between adaptive and
# neutral regions (based on Fst calculations Maria's MA thesis)
############################################################################################################

fst_dir <- "../fst-signals"
fst_files <- list.files(fst_dir, pattern = "*.stats", full.names = TRUE)
fst_data_list <- lapply(fst_files, read.csv)

fst_data_combined <- do.call(rbind, fst_data_list)
fst_data_combined <- fst_data_combined %>%
  dplyr::rename(Chr = scaffold, Start = start, Fst = Fst_S_HUE_T_S_COT_S) %>%
  dplyr::filter(Fst >= 0)

fst_data_combined$chromosome <- dplyr::case_when(
  fst_data_combined$Chr == "CM008230.1_RagTag" ~ "Chr1",
  fst_data_combined$Chr == "CM008233.1_RagTag" ~ "Chr2",
  fst_data_combined$Chr == "CM008234.1_RagTag" ~ "Chr3",
  fst_data_combined$Chr == "CM008235.1_RagTag" ~ "Chr4",
  fst_data_combined$Chr == "CM008236.1_RagTag" ~ "Chr5",
  fst_data_combined$Chr == "CM008237.1_RagTag" ~ "Chr6",
  fst_data_combined$Chr == "CM008238.1_RagTag" ~ "Chr7",
  fst_data_combined$Chr == "CM008239.1_RagTag" ~ "Chr8",
  fst_data_combined$Chr == "CM008240.1_RagTag" ~ "Chr9",
  fst_data_combined$Chr == "CM008231.1_RagTag" ~ "Chr10",
  TRUE ~ NA_character_
  )

fst_data_combined$chromosome <- factor(fst_data_combined$chromosome,
                                       levels = paste0("Chr", 1:10))



# Fst plot trial
library(ggplot2)
ggplot(fst_data_combined, aes(x = Start, y = Fst, color=Fst)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept=0.2) +
  geom_smooth(method = "loess", se = F, span = 0.05, color = "blue") +
  scale_color_gradient(low = "grey", high = "red", guide = "none") +
  facet_wrap(~chromosome, scales = "free_x") + 
  labs(
    title = "Fst Across Chromosomes",
    x = "Gene position (x10^7 bp)",
    y = "Fst"
  ) +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 1), # Rotate and adjust side labels
    legend.position = "bottom"
  ) +
  scale_x_continuous(labels = function(x) x / 1e7)


library(dplyr)
library(purrr)

# Define a helper function to get crossing points per chromosome
find_loess_crossings <- function(df, span = 0.05, threshold = 0.2) {
  df <- df %>% arrange(Start)
  
  loess_fit <- loess(Fst ~ Start, data = df, span = span)
  df$loess_pred <- predict(loess_fit, newdata = df)
  
  # Find where the LOESS prediction crosses the threshold
  crossings <- which(diff(sign(df$loess_pred - threshold)) != 0)
  
  # Interpolate crossing positions
  approx_crossings <- map_dfr(crossings, function(i) {
    x1 <- df$Start[i]
    x2 <- df$Start[i + 1]
    y1 <- df$loess_pred[i]
    y2 <- df$loess_pred[i + 1]
    
    x_cross <- x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
    tibble(Start = x_cross, Fst = threshold)
  })
  
  return(approx_crossings)
}

# Apply per chromosome
crossings_all <- fst_data_combined %>%
  group_by(chromosome) %>%
  group_modify(~ find_loess_crossings(.x)) %>%
  ungroup()
crossings_all
unique(fst_data_combined$chromosome)

crossings_inversions <- crossings_all[-c(4:7, 9, 10, 15, 16, 18, 19, 23:26, 29:32, 34, 35), ]
head(crossings_inversions)
crossings_inversions

ggplot(fst_data_combined, aes(x = Start, y = Fst, color = Fst)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, span = 0.05, color = "#090C02", size = 1) +
  geom_point(data = crossings_inversions, aes(x = Start, y = Fst), color = "#ED6A5E", size = 2) +
  facet_wrap(~chromosome, scales = "free_x") +
  scale_color_gradient(low = "#CCDBDC", high = "#007EA7", guide = "none") +
  scale_x_continuous(labels = function(x) x / 1e7) +
  labs(
    title = "Inversion coordinate filtering: LOESS (span=0.05) and Fst=0.2",
    x = "Position (x10^7 bp)",
    y = "Fst"
  ) +
  theme_minimal()

crossings_inversions$Chr <- dplyr::case_when(
  crossings_inversions$chromosome == "Chr1" ~ "CM008230.1_RagTag",
  crossings_inversions$chromosome == "Chr2" ~ "CM008233.1_RagTag",
  crossings_inversions$chromosome == "Chr3" ~ "CM008234.1_RagTag",
  crossings_inversions$chromosome == "Chr4" ~ "CM008235.1_RagTag",
  crossings_inversions$chromosome == "Chr5" ~ "CM008236.1_RagTag",
  crossings_inversions$chromosome == "Chr6" ~ "CM008237.1_RagTag",
  crossings_inversions$chromosome == "Chr7" ~ "CM008238.1_RagTag",
  crossings_inversions$chromosome == "Chr8" ~ "CM008239.1_RagTag",
  crossings_inversions$chromosome == "Chr9" ~ "CM008240.1_RagTag",
  crossings_inversions$chromosome == "Chr10" ~ "CM008231.1_RagTag",
  TRUE ~ NA_character_
)

# write.table(crossings_inversions, file = "intermediate-output/fst-intervals.tsv", sep = "\t", quote = FALSE, row.names = T)
# write.table(fst_data_combined, file = "intermediate-output/fst-data.tsv", sep = "\t", quote = FALSE, row.names = T)


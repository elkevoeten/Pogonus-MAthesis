library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#####################################################
# small script to get wing size distribution graphs #
#####################################################

winglength <- read.csv("../data/winglength.csv", header = T, dec = ".")
head(winglength)
str(winglength)

winglength <- winglength[-c(1, 2, 8), ]
winglength$ecotype <- as.factor(winglength$ecotype)
winglength$label <- as.factor(winglength$label)
winglength$sex <- as.factor(winglength$sex)
winglength$length_wing <- as.numeric(winglength$length_wing)
winglength$length_elytra <- as.numeric(winglength$length_elytra)

library(ggplot2)
winglength_hist <- ggplot(winglength, aes(x=length_wing, fill = ecotype)) +
  geom_histogram(binwidth = 0.1) +
  scale_fill_manual(values = c("LW" = "#E97132", "SW" = "#4E95D9")) +
  scale_x_continuous(name="Wing length (mm)") +
  scale_y_continuous(name="Number") + 
  labs(fill = "Ecotype") +
  ggtitle("Wing length distribution") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
    )
winglength_hist
ggsave(winglength_hist, 
       filename = "winglength_hist.pdf",
       device = "pdf",
       height = 10, width = 15, units = "cm")

winglength_density <- ggplot(winglength, aes(x = length_wing, color = ecotype)) +
  geom_density(linewidth = 1.2) +  # Creates density lines for each ecotype
  scale_x_continuous(name = "Wing length (mm)") +
  scale_y_continuous(name = "Density") + 
  labs(color = "Ecotype") +  # Custom legend title 
  scale_color_manual(values = c("LW" = "#E97132", "SW" = "#4E95D9")) +
  ggtitle("Wing length distribution") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    legend.position = "bottom"  # Position the legend below
  )
winglength_density

ggsave(winglength_density, 
       filename = "winglength_density.pdf",
       device = "pdf",
       height = 10, width = 15, units = "cm")

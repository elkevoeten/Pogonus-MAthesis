library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#####################################################################################################
# script used to analyse the individual inundation data of the beetles caught in Spain in Sept 2024 #
#####################################################################################################

library(table1)   # For frequency tables
library(car)      # For Type III Anova
library(dplyr)    # For data manipulation
library(e1071)    # For skewness
library(ggplot2)  # For plotting
library(ggpubr)   # For stats on the plots
library(survival) # For survival analysis
library(survminer) # For survival analysis

# DATA LOADING AND EXPLORATION  ---------------------------------------------------------------------------------------------
inundation_sept <- read.csv("../data/ind-inundation-sept.csv", header = T, dec = ".")
head(inundation_sept)
str(inundation_sept)

inundation_sept <- inundation_sept[-c(11, 12, 13, 14, 15, 16), ] #remove failed inundations
inundation_sept$ecotype <- as.factor(inundation_sept$ecotype)
inundation_sept$escape <- ifelse(inundation_sept$time_s == 1200, "no", "yes")
inundation_sept$habitat <- ifelse(inundation_sept$ecotype == "LW", "seasonal", "tidal")
inundation_sept$habitat <- factor(inundation_sept$habitat, levels = c("tidal", "seasonal"))
inundation_sept$ecotype <- as.factor(inundation_sept$ecotype)
inundation_sept$region <- as.factor(inundation_sept$label)
inundation_sept$escape <- as.factor(inundation_sept$escape) # Binomial response (0 = no, 1 = yes)

# frequency table
table1(~ region + habitat | sex, data=inundation_sept)

# visualize the time data
hist(inundation_sept$time_s, main = "Histogram of Time to Resurface", xlab = "Time (seconds)", breaks = 30)

plot(inundation_sept$time_s ~ inundation_sept$ecotype, ylab="Time (s)", xlab="Ecotype",col="lightblue")
stripchart(inundation_sept$time_s ~ inundation_sept$ecotype, vertical = TRUE, method = "jitter", pch = 19, jitter = 0.2, col = "darkblue", add = TRUE)

ggplot(inundation_sept, aes(x = ecotype, y = time_s, color=ecotype)) +
  geom_violin() +
  scale_color_manual(values = c("SW" = "#4E95D9", "LW" = "#E97132")) +
  geom_jitter(width = 0.2, color = "darkblue", size = 1.5) +  # Add jittered points
  labs(x = "Ecotype", y = "Time (s)") +
  theme_minimal()

# MODEL CONSTRUCTION ---------------------------------------------------------------------------------------------
skewness(inundation_sept$time_s)
# Positive skew (> 0) means the data has a right tail, where the majority of the data points are concentrated on the left side
# Okay to proceed with the analysis using models like GLM, which can handle slightly skewed data

# Fit a Gamma GLM with a log link function
# I use a gamma GLM since I have continuous time data (binomial would make sense if i was looking at the binary escape or no escape)
# The log link function is used to ensure that the predicted values stay positive, which is important since the time cannot be negative
glm_model <- glm(time_s ~ ecotype, family = Gamma(link = "log"), data = inundation_sept)
summary(glm_model)

# Perform Type III Anova to test the significance of the factors
# LR, likelihood ratio
# looks at how well adding a specific factor improves the fit of the model
anova_glm <- Anova(glm_model, type = 3, test = "LR")
print(anova_glm)
# ecotype significantly affects the time to resurface -> two ecotypes (SW vs LW) differ in how long it takes them to resurface after being submerged.
# SW beetles took significantly longer to emerge compared to LW beetles

## TRY OUT OTHER MODELS!
# Fit a model with both ecotype and region as predictors
glm_model_region <- glm(time_s ~ ecotype + region, family = Gamma(link = "log"), data = inundation_sept)
AIC(glm_model_region)

# Fit a model with an interaction between ecotype and region
glm_model_interaction <- glm(time_s ~ ecotype * region, family = Gamma(link = "log"), data = inundation_sept)
AIC(glm_model_interaction)

# Fit a Gaussian GLM
glm_model_gaussian <- glm(time_s ~ ecotype, family = gaussian(link = "identity"), data = inundation_sept)
AIC(glm_model_gaussian)

# None of these models are better than the first one (according to their AIC)

# Get the coefficient for 'ecotype' from the GLM model
# A coefficient represents the strength and direction of the effect of that the ecotype has on the outcome (time)
coef_estimate <- coef(glm_model)["ecotypeSW"]

# Exp the coefficient to get the relative effect size (transforms it back to the original scale, since we logged it before)
# the effect size tells me how much the ecotype impacts the response variable (time, in this case)
effect_size <- exp(coef_estimate)
print(effect_size) #2.25, greater than 1 means SW resurfaces later

# Comparing mean values directly, test significance in a simpler way
# Data is not normal, use non-parametric test
wilcox_test_result <- wilcox.test(time_s ~ ecotype, data = inundation_sept)
print(wilcox_test_result)

# SURVIVAL ANALYSIS ----------------------------------------------------------------------------------------------
# Survival analysis is done when the time to event is important, in this case, if beetles have a value of 1200, this means they did not emerge
# If some beetles did not resurface within 20 minutes, their times are "censored"
inundation_sept <- inundation_sept %>%
  mutate(event = ifelse(time_s < 1200, 1, 0))
str(inundation_sept)

# Create a survival object (to use in the model)
surv_obj <- Surv(time = inundation_sept$time_s, event = inundation_sept$event)

# Kaplan-Meier model (NON parametric), estimates the survival probability over time for each ecotype
# It provides an unbiased estimate of the proportion of individuals that remain at risk of resurfacing at each time point
km_fit <- survfit(surv_obj ~ ecotype, data = inundation_sept)

# Plot the survival curve (km plot)
KMsurvivalcurve <- ggsurvplot(km_fit, data = inundation_sept, pval = TRUE,
                              xlab = "Time (seconds)", ylab = "P(staying under water)",
                              #title = "Kaplan-Meier Survival Curve",
                              legend.title = "Habitat", palette = c("#E97132", "#4E95D9"),
                              ggtheme = theme_minimal(base_size = 15) +
                                theme(
                                  legend.title = element_text(face = "bold"),
                                  legend.position = "top",
                                  plot.title = element_text(hjust = 0.5, face = "bold", size = 28),
                                  axis.title.y = element_text(size = 18),
                                  axis.title.x = element_text(size = 18),
                                  axis.text.x = element_text(size = 18),
                                  axis.text.y = element_text(size = 18)
                                ),
                              legend.labs = c("Seasonal", "Tidal"))
KMsurvivalcurve
ggplot_KMsurvivalcurve <- KMsurvivalcurve$plot
ggplot_KMsurvivalcurve <- ggplot_KMsurvivalcurve +
  geom_vline(xintercept = 300, color = "grey", linetype = "solid") +
  guides(fill = guide_legend(title = "Ecotype", title.theme = element_text(face = "bold")))
ggplot_KMsurvivalcurve

ggsave(ggplot_KMsurvivalcurve, 
       filename = "ind-inundation-sept-KMsurvivalcurve.pdf",
       device = "pdf",
       height = 8, width = 10, units = "in")

# Cox proportional hazards model
# tests the significance of ecotype on the hazard rate (rate of emergence over time)
cox_model <- coxph(surv_obj ~ ecotype, data = inundation_sept)
summary(cox_model)

# calculate Hazard Ratio (HR) and confidence intervals for easier interpretation 
# Hazard is the immediate "risk" of coming out of the water for a beetle at a given time
# A higher hazard means a beetle is more likely to resurface quickly
exp(coef(cox_model))  # HR
# ecotypeSW, 0.4182992

# The hazard ratio for ecotypeSW is 0.4183. This value represents the likelihood of resurfacing for SW beetles compared to 
# LW beetles. Since the HR is less than 1, it indicates that SW beetles are less likely to resurface
# The hazard of resurfacing for SW beetles is 58% (1 - 0.4183 x 100) lower than that of LW beetles

# VISUALISATION ---------------------------------------------------------------------------------------------
# boxplot with jitter points
ggplot(inundation_sept, aes(x = ecotype, y = time_s, fill = ecotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot without outliers
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") +  # Add jittered points
  stat_compare_means(label.y = max(inundation_sept$time_s) + 50) +  # Add significance
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Ecotype") +
  labs(title = "Emergence Time by Ecotype") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#E97132", "#4E95D9"))

# violin plot with boxplot on top of it
ggplot(inundation_sept, aes(x = ecotype, y = time_s, fill = ecotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot for density
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +  # Overlay boxplot
  stat_compare_means(label.y = max(inundation_sept$time_s) + 650) +
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Ecotype") +
  labs(title = "Emergence Time by Ecotype") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#E97132", "#4E95D9")) 

# density plot (time is a continuous variable, doesn't need to be shown as steps)
densityplot <- ggplot(inundation_sept, aes(x = time_s, fill = ecotype)) +
  geom_density(alpha = 0.6) +  # Density plot
  scale_x_continuous(name = "Emergence Time (seconds)") +
  scale_y_continuous(name = "Density") +
  labs(title = "Density of Emergence Time by Ecotype") +
  theme_minimal() +
  scale_fill_manual(name = "Ecotype", values = c("#E97132", "#4E95D9")) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
densityplot

ggsave(densityplot, 
       filename = "ind-inundation-sept-density.pdf",
       device = "pdf",
       height = 8, width = 10, units = "in")

# combined (box + violin + jitter)
ind_inundation_september <- ggplot(inundation_sept, aes(x = habitat, y = time_s, fill = habitat)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
  stat_compare_means(label.y = max(inundation_sept$time_s) + 650) +
  scale_y_continuous(
    name = "Emergence Time (seconds)",
    breaks = seq(0, max(inundation_sept$time_s) + 1000, by = 500)
  ) +
  scale_x_discrete(name = "Habitat", labels = c("tidal" = "Tidal", "seasonal" = "Seasonal")) +
  labs(title = "Individual inundations (Spain): September 2024") +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size=12),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.line.y = element_line(color = "darkgrey", size = 0.2)
  )+
  scale_fill_manual(values = c("#4E95D9", "#ffaf7a"))
ind_inundation_september

ggsave(ind_inundation_september, 
       filename = "ind-inundation-sept-violin.pdf",
       device = "pdf",
       height = 8, width = 10, units = "in")







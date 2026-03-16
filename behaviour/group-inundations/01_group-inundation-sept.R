library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#####################################################################################################
# script used to analyse the individual inundation data of the beetles caught in Spain in Sept 2025 #
# Trial 1 and 2: beetles were put back in the wet zone before each flooding event                   #
# Trial 3: beetles were left undisturbed with each flooding event                                   #
#####################################################################################################

# Trial 1 and 2 data (start in wet zone each cycle) --------------------------------------------------------------------------------------------------------
group_inundation_sept <- read.csv("../data/group-inundation-sept-T1-T2.csv", header = T, dec = ".")
head(group_inundation_sept)
str(group_inundation_sept)
group_inundation_sept$ecotype <- as.factor(group_inundation_sept$ecotype)
inundation_T1_T2 <- group_inundation_sept[group_inundation_sept$trial != "T3", ]

#raw data ggplot
library(ggplot2)
inundation_T1_T2_plot <- ggplot(inundation_T1_T2, aes(x = cycle, y = X2H.tot.prop.dry, color = ecotype, group = interaction(label, trial), shape = trial)) + 
  geom_line() +  # Use lines to connect the points
  geom_point(size = 3) +  # Add points for each observation
  scale_color_manual(
    name = "Ecotype", 
    values = c("SW" = "#4E95D9", "LW" = "#E97132"),
    labels = c("LW", "SW")  # Custom labels in legend
  ) +  
  scale_shape_manual(
    name = "Trial", 
    values = c(16, 17),  # Adjust based on trial levels
    labels = c("T1", "T2")  # Custom labels in legend
  ) +  
  labs(x = "Day", y = "Proportion on dry part after 2h of inundation", title = "Change in reaction to flooding") + 
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(inundation_T1_T2_plot, 
       filename = "group-inundation-sept-T1-T2.pdf",
       device = "pdf")

#Statistics
#Exploratory data analysis:
summary(inundation_T1_T2)
ggplot(inundation_T1_T2, aes(x = ecotype, y = X2H.tot.prop.dry, fill = ecotype)) +
  geom_boxplot() +
  facet_wrap(~ trial) +
  theme_minimal()
ggplot(inundation_T1_T2, aes(x = cycle, y = X2H.tot.prop.dry, color = ecotype, group = ecotype)) +
  geom_line(stat = "summary", fun = mean) +
  geom_point(stat = "summary", fun = mean, size = 3) +
  facet_wrap(~ trial) +
  theme_minimal()

#Statistical model
library(lme4)
model_simple <- lm(X2H.tot.prop.dry ~ cycle + ecotype, data = inundation_T1_T2)
summary(model_simple)
# This p-value when comparing habitat type (SW vs LW) is significant (p < 0.01), indicating that the difference between the 
# SW and LW habitats is highly significant (p = 0.0029).

model_interaction <- lm(X2H.tot.prop.dry ~ cycle * ecotype, data = inundation_T1_T2)
summary(model_interaction)
# the interaction between time and habitat is not significant (p = 0.58695)
# This means that the interaction between time and habitat has no effect 
# on the proportion of individuals on the dry side. 

model_int_fixed <- lmer(X2H.tot.prop.dry ~ ecotype * cycle + (1 | label), data = inundation_T1_T2)
summary(model_int_fixed)

AIC(model_simple)       # AIC for the simpler model, -2.898342 -> better model
AIC(model_interaction)  # AIC for the model with interaction, -1.299373  
AIC(model_int_fixed)    # AIC for model with interaction & fixed effect, 16.372

#assumptions
plot(residuals(model_simple))
hist(residuals(model_simple))
qqnorm(residuals(model_simple))
qqline(residuals(model_simple), col = "red")
shapiro.test(residuals(model_simple)) #W = 0.94231, p-value = 0.3466


# Create predictions from the model
inundation_T1_T2$predicted <- predict(model_simple)
# Plot observed vs predicted values
ggplot(inundation_T1_T2, aes(x = cycle, y = X2H.tot.prop.dry, color = ecotype)) +
  geom_point() +
  geom_line(aes(y = predicted), size = 1, linetype = "dashed") +
  labs(
    title = "Observed vs predicted proportions (lm(X2H.tot.prop.dry ~ cycle + ecotype, data = inundation_T1_T2))",
    x = "Time (days)",
    y = "Proportion on Dry Side",
    color = "Ecotype"
  ) +
  theme_minimal()

ggplot(inundation_T1_T2, aes(x = X2H.tot.prop.dry)) + 
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  theme_minimal()
#A histogram helps you visually inspect the distribution of the data. If the data looks approximately 
#bell-shaped and symmetric, it suggests normality (parametric).
#If the histogram is skewed or shows multiple peaks, the data might not follow a normal distribution.
qqnorm(inundation_T1_T2$X2H.tot.prop.dry)
qqline(inundation_T1_T2$X2H.tot.prop.dry, col = "red")
#If the points lie along the red line, the data likely follows a normal distribution.
#If the points deviate significantly from the line (e.g., curves or outliers), the data may not be normal.
shapiro.test(inundation_T1_T2$X2H.tot.prop.dry)
#Since this p-value is greater than 0.05, you fail to reject the null hypothesis, 
#suggesting that the data is consistent with a normal distribution (i.e., parametric).
#IF THE DATA DOES NOT FOLLOW A NORMAL DISTRIBUTION, USE NON-PARAMETRIC TESTS

# Perform ANOVA to test for differences between ecotypes and cycles
anova_results <- aov(X2H.tot.prop.dry ~ ecotype + cycle, data = inundation_T1_T2)
# View the results of the ANOVA
summary(anova_results)
summary(model_simple)

# assumptions 
# Plot residuals for linear regression model
par(mfrow = c(2, 2))  # Set up for multiple plots
plot(model_simple)  # This will show diagnostic plots

# You can also check for homogeneity of variance using the Breusch-Pagan test:
install.packages("lmtest")
library(lmtest)
bptest(model_simple)

# visualisation
library(effects)
plot(predictorEffects(model_simple, ~ ecotype + cycle), multiline = TRUE)
plot(allEffects(model_simple), multiline = TRUE, confint = list(style = "auto"))

# Load required packages
library(ggplot2)
library(dplyr)

# Generate model predictions with confidence intervals
predictions <- predict(model_simple, newdata = inundation_T1_T2, interval = "confidence")

# Add predictions and confidence intervals to the dataset
inundation_T1_T2 <- inundation_T1_T2 %>%
  mutate(predicted = predictions[, "fit"],
         lower = predictions[, "lwr"],
         upper = predictions[, "upr"])

# Now plot the data with raw points, predicted values, and confidence intervals
inundation_T1_T2_plot2 <- ggplot(inundation_T1_T2, aes(x = cycle, y = X2H.tot.prop.dry, color = ecotype)) +
  geom_point(size = 3, alpha = 0.6, aes(shape = trial)) +  
  geom_line(aes(x = cycle, y = X2H.tot.prop.dry, color = ecotype, group = interaction(ecotype, trial)), 
            size = 0.5, alpha = 0.5) +
  geom_line(aes(x = cycle, y = predicted), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ecotype), alpha = 0.2, color = NA) +
  scale_color_manual(name = "Ecotype", values = c("SW" = "#4E95D9", "LW" = "#E97132")) +  
  scale_fill_manual(name = "Ecotype", values = c("SW" = "#4E95D9", "LW" = "#E97132")) +  
  scale_shape_manual(
    name = "Trial", 
    values = c(16, 17),  
    labels = c("T1", "T2"))+
  labs(x = "Cycle", y = "Proportion on dry area after 2h", title = "Prediction of the change in initial response to inundation") +
  annotate("text", x = 3, y = 0.85, label = "ANOVA, ecotype: \np = 0.003", size = 4, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(inundation_T1_T2_plot2, 
       filename = "group-inundation-sept-T1-T2-pred.pdf",
       device = "pdf",
       height = 6, width = 8, bg = "white")

# Trial 3 data (no disturbance) ----------------------------------------------------------------------------------------------------------
# data loading and exploration
inundation_T3 <- read.csv("../data/group-inundation-sept-T3.csv", header=T)
head(inundation_T3)
str(inundation_T3)
inundation_T3$ecotype <- as.factor(inundation_T3$ecotype)

inundation_T3_raw <- ggplot(inundation_T3, aes(x=time_h, y=prop.dry, color=ecotype)) +
  annotate("rect", xmin = 0, xmax = 2, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 24, xmax = 26, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 48, xmax = 50, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 72, xmax = 74, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 96, xmax = 98, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  geom_line() +
  geom_point(size = 1) +
  scale_color_manual(name = "Ecotype", values = c("SW" = "#4E95D9", "LW" = "#E97132")) +
  labs(x = "Time (hours)", y = "Proportion of beetles on dry part", title = "Large scale behaviour without disturbance") +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
    )
inundation_T3_raw

ggsave(inundation_T3_raw, 
       filename = "group-inundation-sept-T3-raw.pdf",
       device = "pdf",
       height = 5, width = 4, bg = "white")

# T3 statistics ----------------------------------------------------------------------------------------------------------
df <- as.data.frame(inundation_T3)
colnames(df)[2] <- "habitat"
df$time_h <- as.numeric(df$time_h)

# Simple model without interaction or random effects
model_simple <- lm(prop.dry ~ time_h + habitat, data = df)
summary(model_simple)
# This p-value when comparing habitat type (SW vs LW) is significant (p < 0.01), indicating that the difference between the 
# SW and LW habitats is highly significant (p = 0.002897).

# Model without random effects + interaction between time and habitat
model_interaction <- lm(prop.dry ~ time_h * habitat, data = df)
summary(model_interaction)
# the interaction between time and habitat is highly significant (p = 0.00576)
# This means that the interaction between time and habitat has a significant effect 
# on the proportion of individuals on the dry side. 

# Calculate AIC for each model
AIC(model_simple)     # AIC for the simpler model, 9.171824
AIC(model_interaction)  # AIC for the model with interaction, 1.351694 -> better model 

plot(residuals(model_interaction))
hist(residuals(model_interaction))
qqnorm(residuals(model_interaction))
qqline(residuals(model_interaction), col = "red")
# model seems to fit the data quite well

# Create predictions from the model
df$predicted <- predict(model_interaction)
# Plot observed vs predicted values
ggplot(df, aes(x = time_h, y = prop.dry, color = habitat)) +
  geom_point() +
  geom_line(aes(y = predicted), size = 1, linetype = "dashed") +
  labs(
    title = "Observed vs predicted proportions (lm(prop.dry ~ time_h * habitat)",
    x = "Time (hours)",
    y = "Proportion on Dry Side",
    color = "Habitat"
  ) +
  theme_minimal()

# Line plot
ggplot(df, aes(x = time_h, y = prop.dry, group = interaction(habitat, label), color = habitat)) +
  geom_vline(xintercept = c(0, 24, 48, 72, 96), color = "lightblue", size = 3, alpha = 0.3) +
  geom_vline(xintercept = c(2, 26, 50, 74, 98), color = "lightcoral", size = 1, alpha = 0.2) +
  # Main plot elements
  geom_line(size = 1) +
  geom_point(size = 2) +
  # Labels and theme
  labs(
    title = "Spatial sorting of ecotypes in response to flooding",
    x = "Time (hours)",
    y = "Proportion on Dry Side",
    color = "Habitat",
    linetype = "Label"
  ) +
  theme_minimal()

group_inundation_T3_plot <- ggplot(df, aes(x=time_h, y=prop.dry, group = interaction(habitat, label), color=habitat)) +
  annotate("rect", xmin = 0, xmax = 2, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 24, xmax = 26, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 48, xmax = 50, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 72, xmax = 74, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  annotate("rect", xmin = 96, xmax = 98, ymin = -Inf, ymax = Inf, fill = "#DCEAF7") +
  geom_vline(xintercept = c(2, 26, 50, 74, 98), color = "lightcoral", size = 1, alpha = 0.6, linetype = "solid") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("LW" = "#E97132", "SW" = "#4E95D9")
  ) +
  guides(color = guide_legend(title = "Ecotype"))+
  labs(
    x = "Time (hours)", 
    y = "Proportion of beetles in dry zone", 
    title = "Spatial sorting of ecotypes in response to flooding",
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size=15),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )
group_inundation_T3_plot

ggsave(group_inundation_T3_plot, 
       filename = "group-inundation-sept-T3-plot.pdf",
       device = "pdf",
       height = 6.5, width = 10, bg = "white")


# I don't think the data is parametric, let's check!
ggplot(df, aes(x = prop.dry)) + 
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  theme_minimal()
#A histogram helps you visually inspect the distribution of the data. If the data looks approximately 
#bell-shaped and symmetric, it suggests normality (parametric).
#If the histogram is skewed or shows multiple peaks, the data might not follow a normal distribution.
qqnorm(df$prop.dry)
qqline(df$prop.dry, col = "red")
#If the points lie along the red line, the data likely follows a normal distribution.
#If the points deviate significantly from the line (e.g., curves or outliers), the data may not be normal.
shapiro.test(df$prop.dry)
#Since this p-value is greater than 0.05, you fail to reject the null hypothesis, 
#suggesting that the data is consistent with a normal distribution (i.e., parametric).
# THE DATA DOES NOT FOLLOW A NORMAL DISTRIBUTION, USE NON-PARAMETRIC TESTS

# Perform Wilcoxon rank-sum test (Mann-Whitney U test) since we just determined that the data is non-parametric
wilcox_test_result <- wilcox.test(prop.dry ~ habitat, data = df)
p_value <- wilcox_test_result$p.value
p_value

# Line plot with stats
ggplot(df, aes(x = time_h, y = prop.dry, group = interaction(habitat, label), color = habitat)) +
  geom_vline(xintercept = c(0, 24, 48, 72, 96), color = "lightblue", size = 3, alpha = 0.3) +
  geom_vline(xintercept = c(2, 26, 50, 74, 98), color = "lightcoral", size = 1, alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Spatial sorting of ecotypes in response to flooding",
    x = "Time (hours)",
    y = "Proportion of individuals on dry side",
    color = "Habitat",
    linetype = "Habitat"
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = 14,  # Adjust x-position as needed
    y = 1,   # Adjust y-position as needed
    label = paste("p-value =", signif(p_value, 3)),
    color = "black",
    size = 3
  )


# Create a table of counts for each habitat and each condition (dry vs. wet)
table_data <- table(df$habitat, df$prop.dry > 0.5)
# Proportions test
prop_test <- prop.test(table_data)
print(prop_test)
# statistically significant difference in the proportions between the two groups (p-value = 0.02248)

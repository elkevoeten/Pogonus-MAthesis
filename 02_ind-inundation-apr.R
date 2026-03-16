library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

####################################################################################################
# script used to analyse the individual inundation data of the beetles caught in Spain in Apr 2025 #
####################################################################################################

library(table1) # For frequency tables
library(car) # For Type III Anova
library(dplyr) # For data manipulation
library(e1071) # For skewness
library(ggplot2) # For plotting
library(ggpubr) # For stats on the plots
library(survival) # For survival analysis
library(survminer) # For survival analysis

# DATA LOADING AND EXPLORATION  ---------------------------------------------------------------------------------------------
inundation_april <- read.csv("../data/ind-inundation-apr.csv", header = T, dec = ".")
head(inundation_april)
str(inundation_april)
nrow(inundation_april)
inundation_april$species <- as.factor(inundation_april$species)
inundation_april$habitat <- as.factor(inundation_april$habitat)
inundation_april$sex <- as.factor(inundation_april$sex)
inundation_april$label <- as.factor(inundation_april$label)
inundation_april$species_ET <- ifelse(
  inundation_april$species == "chalceus",
  ifelse(inundation_april$habitat == "seasonal", "P. chalceus (LW)", "P. chalceus (SW)"),
  paste0("P. ", inundation_april$species)
)
inundation_april$species_ET <- as.factor(inundation_april$species_ET)

table1(~ label + habitat | sex, data=inundation_april)

plot(inundation_april$time_s ~ inundation_april$label, ylab="Time (s)", xlab="Label",col="lightblue")
stripchart(inundation_april$time_s ~ inundation_april$habitat, vertical = TRUE, method = "jitter", pch = 19, jitter = 0.2, col = "darkblue", add = TRUE)

ggplot(inundation_april, aes(x = label, y = time_s, color=label)) +
  geom_violin() +
  #scale_color_manual(values = c("tidal" = "#4E95D9", "seasonal" = "#E97132")) +
  geom_jitter(width = 0.2, color = "darkblue", size = 1.5) +  # Add jittered points
  labs(x = "Habitat", y = "Time (s)") +
  theme_minimal()

ind_inundation_apr <- ggplot(inundation_april, aes(x = habitat, y = time_s, fill = habitat)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, color = "black") +
  stat_compare_means(label.y = max(inundation_april$time_s) + 100) +
  scale_y_continuous(name = "Emergence Time (seconds)") +
  scale_x_discrete(name = "Habitat") +
  labs(title = "Individual inundations (Spain, April 2025)") +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size=28),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  )+
  scale_fill_manual(values = c("#E97132", "#4E95D9"))
ind_inundation_apr

ggsave(ind_inundation_apr, 
       filename = "ind-inundation-apr-violin.pdf",
       device = "pdf",
       height = 8, width = 10, units = "in")

levels(inundation_april$label)
inundation_april$label <- factor(inundation_april$label, levels = c("BAR", "TROC", "CET", "GILV", "LITT", "MER", "SMAR"))

levels(inundation_april$species_ET)
inundation_april$species_ET <- factor(inundation_april$species_ET, levels = c("P. chalceus (SW)", "P. chalceus (LW)", "P. gilvipes", "P. littoralis", "P. meridionalis", "P. smaragdinus"))

ind_inundation_apr_species <- ggplot(inundation_april, aes(x = species_ET, y = time_s, fill = species_ET)) +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "darkgrey") +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.8, color = "black") +
  #stat_compare_means(label.y = max(inundation_april$time_s) + 100) +
  scale_y_continuous(name = "Emergence time (seconds)") +
  scale_x_discrete(name = "Species (ecotype)") +
  labs(title = "April 2025") +
  theme_minimal() +
  scale_fill_manual(values = c("#4E95D9", "#E97132","#ffaf7a", "#ffaf7a" ,"#ffaf7a", "#ffaf7a")) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size=12),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "italic", size = 8),
    axis.text.y = element_text(size = 8),
    axis.line.y = element_line(color = "darkgrey", size = 0.2)
  )
ind_inundation_apr_species
ggsave(ind_inundation_apr_species, 
       filename = "ind-inundation-apr-species.pdf",
       device = "pdf",
       height = 12, width = 20, units = "cm")

kruskal.test(time_s ~ label, data = inundation_april)

# MODEL AND ASSUMPTIONS ---------------------------------------------------------------------------------------------
# Visualize the time data
hist(inundation_april$time_s, main = "Histogram of Time to Resurface", xlab = "Time (seconds)", breaks = 30)

# Remove the rows that contain the NA
inundation_april <- na.omit(inundation_april)

# Check for skewness
skewness(inundation_april$time_s)
# Positivinundation# Positive skew (> 0) means the data has a right tail, where the majority of the data points are concentrated on the left side
# Okay to proceed with the analysis using models like GLM, which can handle slightly skewed data

# Fit a Gamma GLM with a log link function
# I use a gamma GLM since I have continuous time data (binomial would make sense if i was looking at the binary escape or no escape)
# The log link function is used to ensure that the predicted values stay positive, which is important since the time cannot be negative
glm_model <- glm(time_s ~ habitat, family = Gamma(link = "log"), data = inundation_april)

# Summary of the model
summary(glm_model) #AIC 1102.5

# Perform Type III Anova to test the significance of the factors
# LR, likelihood ratio
# looks at how well adding a specific factor improves the fit of the model
anova_glm <- Anova(glm_model, type = 3, test = "LR")
print(anova_glm)
#  ecotype significantly affects the time to resurface -> two ecotypes (SW vs LW) differ in how long it takes them to resurface after being submerged.
# SW beetles took significantly longer to emerge compared to LW beetles

## TRY OUT OTHER MODELS!
# Fit a model with both ecotype and region as predictors
glm_model_region <- glm(time_s ~ habitat + sex, family = Gamma(link = "log"), data = inundation_april)
AIC(glm_model_region) #1103.883

# Fit a model with an interaction between ecotype and region
glm_model_interaction <- glm(time_s ~ habitat + label, family = Gamma(link = "log"), data = inundation_april)
AIC(glm_model_interaction)

# Fit a Gaussian GLM
glm_model_gaussian <- glm(time_s ~ habitat, family = gaussian(link = "identity"), data = inundation_april)
AIC(glm_model_gaussian)

# None of these models are better than the first one (according to their AIC)

# Get the coefficient for 'ecotype' from the GLM model
# A coefficient represents the strength and direction of the effect of that the ecotype has on the outcome (time)
coef_estimate <- coef(glm_model)["habitattidal"]

# Exp the coefficient to get the relative effect size (transforms it back to the original scale, since we logged it before)
# the effect size tells me how much the ecotype impacts the response variable (time, in this case)
effect_size <- exp(coef_estimate)
print(effect_size) #1.776452, greater than 1 means tidal resurfaces later

# Comparing mean values directly, test significance in a simpler way
# Data is not normal, use non-parametric test
wilcox_test_result <- wilcox.test(time_s ~ habitat, data = inundation_april)
print(wilcox_test_result)

# SURVIVAL ANALYSIS ----------------------------------------------------------------------------------------------
library(survival)
library(flexsurv)
inundation_april$status <- ifelse(inundation_april$time_s == 1200, 0, 1) #setting the status to 0 if there was no escape at 1200s, otherwise 1 for escape
cox_fit <- coxph(Surv(time_s, status) ~ habitat + sex + label, data = inundation_april)
summary(cox_fit)
cox_zph <- cox.zph(cox_fit)
print(cox_zph) #no p<0.05, so no violation of the proportional hazards assumption
plot(cox_zph)

fit <- survreg(Surv(time_s, status, type = "right") ~ habitat, data=inundation_april)

fit_exp_flexsurv <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat, dist = "exponential", data = inundation_april)
fit_wei_flexsurv <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat, dist = "weibull", data = inundation_april)
fit_gen_flexsurv <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat, dist = "gengamma", data = inundation_april)
fit_logn_flexsurv <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat, dist = "lnorm", data = inundation_april)
fit_logl_flexsurv <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat, dist = "llogis", data = inundation_april)

aicc_values <- data.frame(
  Model = c("Exponential", "Weibull", "Generalized Gamma", "Log Normal", "Log Logistic"),
  AICc = c(AICc(fit_exp_flexsurv), AICc(fit_wei_flexsurv), AICc(fit_gen_flexsurv), 
           AICc(fit_logn_flexsurv), AICc(fit_logl_flexsurv))
)
print(aicc_values)
#               Model     AICc
# 1       Exponential 2425.326
# 2           Weibull 2388.124
# 3 Generalized Gamma 2275.101 <<<
# 4        Log Normal 2316.640
# 5      Log Logistic 2321.826

fit_gen_flexsurv2 <- flexsurvreg(Surv(time_s, status, type = "right") ~ label, dist = "gengamma", data = inundation_april)
fit_gen_flexsurv3 <- flexsurvreg(Surv(time_s, status, type = "right") ~ habitat + sex, dist = "gengamma", data = inundation_april)
fit_gen_flexsurv4 <- flexsurvreg(Surv(time_s, status, type = "right") ~ label + sex, dist = "gengamma", data = inundation_april)

AICc(fit_gen_flexsurv2) # 2275.478
AICc(fit_gen_flexsurv3) # 2274.311
AICc(fit_gen_flexsurv4) # 2271.375 <<<<<

#### Plot mean non-escape of beetles and hazard in function of time ####
plot(fit_gen_flexsurv4, ylab = "Not escaped beetles", xlab = "time (s)", main = "gengamma")
plot(fit_gen_flexsurv4, ylab = "Hazard", xlab = "time (s)", main = "gengamma", type="hazard")

#### Inspect summary coefficient table & confidence intervals ####
# Which of the covariates (habitat, label) significantly affects escape when used in a predictive model?
fit_gen_flexsurv4$res 
fit_gen_flexsurv$res

#### Check the quality of the fit based on the deviance #### --> couldn't use deviance, so did coxsnell and there is a pattern, suggesting the fit is not ok???
cox = residuals(fit_gen_flexsurv4,type="coxsnell")
par(mfrow = c(1, 2))
for (f in c("time_s","label")) {
  plot(inundation_april[,f],cox,xlab=f,ylab="coxsnell") 
  abline(h = 0, lty = 2)
}
par(mfrow = c(1, 1))







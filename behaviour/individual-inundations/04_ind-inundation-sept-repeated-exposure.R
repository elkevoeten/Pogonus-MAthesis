library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#####################################################################################################
# script used to analyse the individual inundation data of repeated exposure of the beetles caught  # 
# in Spain in Sept 2024                                                                             #
#####################################################################################################

library(readxl)
repeated_exposure <- read_xlsx(path = "../data/ind-inundation-repeats.xlsx")

library(ggplot2)
ggplot(repeated_exposure, aes(x = trial, y = time_seconds, group = individual_id, color = ecotype)) +
  geom_line() +
  geom_point() +
  labs(title = 'Change in Time Over Multiple Trials',
       x = 'Trial',
       y = 'Time (seconds)') +
  theme_minimal()

library(lme4)
fit <- lmer(time_seconds ~ ecotype * trial + (1 | individual_id), repeated_exposure)
fit2 <- lmer(time_seconds ~ ecotype * trial + (1 | individual_id) + (1 | trial), repeated_exposure)

library(effects)
plot(allEffects(fit, residuals= T))
plot(allEffects(fit), multiline=T)
effects_data <- as.data.frame(allEffects(fit))

library(stats)
anova(fit)
summary(fit)

library(effects)
plot(allEffects(fit))

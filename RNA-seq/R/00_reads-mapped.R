# Script to analyse the amount of mapped reads and whether it differs significantly between ecotypes
# -------------------------------------------------------------------------------------------------------------

# Load required packages
if (!require("car")) install.packages("car", dependencies = TRUE)
if (!require("ggpubr")) install.packages("ggpubr", dependencies = TRUE)

library(car)
library(ggpubr)

# Step 1: Create data frame
mapping <- data.frame(
  sample = c("BAR2-01", "BAR2-02", "BAR2-03", "BAR2-04", "BAR2-05", "BAR2-06",
             "BAR4-01", "BAR4-02", "BAR4-03", "BAR4-04", "BAR4-05", "BAR4-06"),
  reads_mapped = c(91080754, 110697248, 100546210, 107391802, 109538760, 110404836,
                   71058604, 65377320, 67358594, 70423050, 72864122, 72016274),
  clean_reads = c(49851927, 60053796, 55628859, 60178967, 60094900, 60102911,
                  60118153, 54941753, 60043103, 60017527, 60055642, 60105182),
  clean_base = c(14955578100, 18016138800, 16688657700, 18053690100, 18028470000, 18030873300,
                 18035445900, 16482525900, 18012930900, 18005258100, 18016692600, 18031554600)
)

# Step 2: Split groups
bar2 <- mapping$reads_mapped[grep("^BAR2", mapping$sample)]
bar4 <- mapping$reads_mapped[grep("^BAR4", mapping$sample)]

# Combine into one data frame for assumption testing
group <- factor(c(rep("BAR2", length(bar2)), rep("BAR4", length(bar4))))
reads_mapped <- c(bar2, bar4)
data <- data.frame(group, reads_mapped)

# Step 3: Test normality (Shapiro-Wilk)
shapiro_bar2 <- shapiro.test(bar2)
shapiro_bar4 <- shapiro.test(bar4)

# Step 4: Test homogeneity of variance (Levene's Test)
levene <- leveneTest(reads_mapped ~ group, data = data)

# Step 5: Q-Q plots to visually assess normality
par(mfrow = c(1, 2))
qqnorm(bar2, main = "Q-Q Plot for BAR2")
qqline(bar2)
qqnorm(bar4, main = "Q-Q Plot for BAR4")
qqline(bar4)

# Step 6: Perform Welch’s t-test
t_test_result <- t.test(bar2, bar4)

# Step 7: Output results
cat("Shapiro-Wilk Test for BAR2:\n")
print(shapiro_bar2)

cat("\nShapiro-Wilk Test for BAR4:\n")
print(shapiro_bar4)

cat("\nLevene's Test for Homogeneity of Variance:\n")
print(levene)

cat("\nWelch's t-test Result:\n")
print(t_test_result)

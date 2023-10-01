library(lme4)
library(ggplot2)

df <- read.csv("parsed_data_all.csv")

# Fit Linear Mixed Effects Model
fit_lmm <- lmer(Volume ~ Time + (1 | PatientID), data = df)

summary(fit_lmm)

#--------------------------

# Predicted values
df$Predicted <- predict(fit_lmm, re.form = NA)

# Bootstrapped confidence intervals
boot_fit <- bootMer(fit_lmm, FUN = function(x) fixef(x), nsim = 1000)
pred_ci <- apply(boot_fit$t, 2, function(x) quantile(x, c(0.025, 0.975)))

ggplot(df, aes(x = Time, y = Volume)) +
     geom_line(aes(color = PatientID, group = PatientID)) + # Individual growth curves
     geom_point(aes(color = PatientID)) + # Data points
     geom_line(aes(y = Predicted), size = 1.5) + # LMM Prediction
     geom_ribbon(aes(ymin = pred_ci[1], ymax = pred_ci[2]), alpha = 0.5) + # Confidence interval
     labs(
          title = "Tumor Growth Over Time",
          x = "Time",
          y = "Tumor Volume",
          color = "Patient ID"
     ) +
     theme_minimal()







library(lme4)
library(ggplot2)
library(dplyr)

# Read the data
df <- read.csv('parsed_data_all.csv')

# Fit the Linear Mixed Effects Model
fit_lmm <- lmer(Volume ~ Time + (1|PatientID), data = df)

# Get unique time points
unique_times <- unique(df$Time)

# Initialize an empty data frame to store CI and predicted values
bootstrap_ci <- data.frame()

# Perform bootstrapping for each unique time point
for(t in unique_times) {
  temp_df <- df[df$Time == t, ]
  
  # Skip iteration if only one observation
  if(nrow(temp_df) <= 1) next
  
  boot_fit <- bootMer(fit_lmm, FUN=function(x) predict(x, newdata = temp_df, re.form = NA), nsim=100)
  
  # Get CI
  ci_low <- quantile(boot_fit$t, 0.025)
  ci_high <- quantile(boot_fit$t, 0.975)
  predicted <- mean(boot_fit$t)
  
  # Append to final data frame
  bootstrap_ci <- rbind(bootstrap_ci, data.frame(Time = t, Predicted = predicted, CI_low = ci_low, CI_high = ci_high))
}

# Plot
ggplot(df, aes(x=Time, y=Volume)) +
  geom_line(aes(group=PatientID, color=PatientID), alpha=0.5) + # Individual growth curves
  geom_point(aes(color=PatientID), alpha=0.5) + # Data points
  geom_line(data=bootstrap_ci, aes(y=Predicted), color="black", size=1.5) + # LMM Prediction
  geom_ribbon(data=bootstrap_ci, aes(y=Predicted, ymin=CI_low, ymax=CI_high), fill="grey", alpha=0.5) + # Confidence interval
  labs(title="Tumor Growth Over Time",
       x="Time",
       y="Tumor Volume",
       color="Patient ID") +
  theme_minimal()

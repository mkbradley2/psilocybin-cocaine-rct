# 03. Survival Analysis.r
# Psilocybin for Cocaine Use Disorder RCT - Survival Analysis and Post-hoc Power Simulations
# Author: Melissa Bradley
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-04-05
#
# Purpose:
#   Conduct survival analysis on time to lapse using Cox models and Kaplan-Meier curves.
#   Assess statistical power post-hoc through simulation, including Firth correction 
#   for small-sample bias (correction not included in final manuscript, see below).

library(tidyverse)   # Wrangling
library(survival)    # Survival models
library(survminer)   # KM plots
library(broom)       # Model tidying
library(coxphf)      # Firth’s correction


# Load data
cocaine_data <- read_sav("Peters data testing SPSS.sav")
time_vars <- read_sav("Observation Days.sav")
status_and_days <- read_sav("Otto survival.sav")

# ================================================
# Data Cleaning for Modeling
# ================================================

cocaine_data_filtered <- cocaine_data %>%
  left_join(select(status_and_days, ID_Sub, STATUS, DAYS), by = "ID_Sub") %>%
  filter(ID_Sub %in% time_vars$ID_Sub, !is.na(Time_Recoded)) %>%
  mutate(
    status_recoded = case_when(
      STATUS == 0 ~ 1,
      STATUS == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    condition_recoded = case_when(
      Condition == 2 ~ 1,
      Condition == 1 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE)

specific_columns <- time_vars %>%
  group_by(ID_Sub) %>%
  mutate(days_under_obs = sum(Prescreen, Preparation, Drug, Integration, na.rm = TRUE)) %>%
  select(ID_Sub, Prescreen, Preparation, Drug, Integration, days_under_obs) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup()

combined_data <- cocaine_data_filtered %>%
  left_join(specific_columns, by = "ID_Sub") %>%
  distinct(ID_Sub, .keep_all = TRUE)

cocaine_data_survival <- combined_data %>%
  mutate(
    days_recoded = case_when(
      STATUS == 1 ~ Integration,
      ID_Drug == 1027 ~ Integration + 90 + 1,
      ID_Drug == 1030 ~ Integration + 90 + 1,
      ID_Drug == 1037 ~ Integration + 1,
      ID_Drug == 1012 ~ 0,
      TRUE ~ if_else(!is.na(DAYS), DAYS + 1, NA_real_)
    )
  ) %>%
  select(ID_Sub, Gender, Ethnicity, Race, Condition, days_recoded, condition_recoded, status_recoded, STATUS, Integration)

# ================================================
# Survival Analysis
# ================================================

km_fit <- survfit(Surv(time = cocaine_data_survival$days_recoded, event = cocaine_data_survival$status_recoded) ~ condition_recoded, data = cocaine_data_survival)
logrank_test <- survdiff(Surv(days_recoded, status_recoded) ~ condition_recoded, data = cocaine_data_survival)
chisq_stat <- logrank_test$chisq
df <- length(logrank_test$n) - 1
p_val <- 1 - pchisq(chisq_stat, df)

cat("Log-Rank Test:\n")
cat("Chi-square =", round(chisq_stat, 3), "\n")
cat("df =", df, "\n")
cat("p-value =", formatC(p_val, digits = 6, format = "f"), "\n")

summary(km_fit)

survival_plot <- ggsurvplot(
  km_fit, 
  data = cocaine_data_survival, 
  pval = TRUE,
  pval.coord = c(0, 0.03),
  conf.int = TRUE,
  xlab = "Days",
  ylab = "Proportion Without Lapse (%)",
  surv.scale = "percent",
  legend.title = "",
  legend.labs = c("Placebo", "Psilocybin"),
  legend = "right",
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "right",
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.8),
      axis.line.y = element_line(color = "black", linewidth = 0.8)
    )
)

ggsave(
  filename = "cocaine_survival_curve_notitle.png",
  plot = survival_plot$plot,
  width = 14,
  height = 7.5,
  dpi = 300
)

surv_object <- Surv(time = cocaine_data_survival$days_recoded, event = cocaine_data_survival$status_recoded)

cox_model_unadjusted <- coxph(surv_object ~ condition_recoded, data = cocaine_data_survival)
cox_model_adjusted <- coxph(surv_object ~ condition_recoded + Gender, data = cocaine_data_survival)

cox_summary_unadjusted <- tidy(cox_model_unadjusted) %>% mutate(Model = "Unadjusted")
cox_summary_adjusted <- tidy(cox_model_adjusted) %>% mutate(Model = "Adjusted for Gender")

cox_summaries <- bind_rows(cox_summary_unadjusted, cox_summary_adjusted) %>%
  select(Model, term, estimate, std.error, statistic, p.value) %>%
  rename(Term = term, Estimate = estimate, `Std. Error` = std.error, `Z-Statistic` = statistic, `P-Value` = p.value) %>%
  mutate(
    Estimate = round(Estimate, 3),
    `Std. Error` = round(`Std. Error`, 3),
    `Z-Statistic` = round(`Z-Statistic`, 3),
    `P-Value` = round(`P-Value`, 3),
    `Hazard Ratio` = round(exp(Estimate), 3),
    `HR Lower CI` = round(exp(Estimate - 1.96 * `Std. Error`), 3),
    `HR Upper CI` = round(exp(Estimate + 1.96 * `Std. Error`), 3)
  )

print(cox_summaries)

# Proportion of participants who had an event (lapse)
event_rate <- mean(cocaine_data_survival$status_recoded)

# Observed hazard ratio from Cox model
obs_hr <- exp(coef(cox_model_unadjusted)["condition_recoded"])

simulate_cox_data <- function(n, hr, event_rate, max_followup = 180) {
  treatment <- rbinom(n, 1, 0.5)
  lambda0 <- -log(1 - event_rate) / max_followup
  hazard <- lambda0 * exp(log(hr) * treatment)
  time <- rexp(n, hazard)
  time <- pmin(time, max_followup)
  status <- as.integer(time < max_followup)
  data.frame(time = time, status = status, treatment = treatment)
}

set.seed(2025)
n_sim <- 1000
n <- 40
p_vals <- replicate(n_sim, {
  sim_data <- simulate_cox_data(n = n, hr = obs_hr, event_rate = event_rate)
  sim_model <- coxph(Surv(time, status) ~ treatment, data = sim_data)
  summary(sim_model)$coefficients["treatment", "Pr(>|z|)"]
})

simulated_power <- mean(p_vals < 0.05)
simulated_power

# Cox Regression w/ Firth's Correction for small samples:
# Although Firth’s correction is often recommended in small samples to reduce 
# bias in hazard ratio estimation, it can substantially increase variance and 
# reduce statistical power due to conservative shrinkage. In this study, the 
# Firth-corrected model yielded even lower simulated power than the uncorrected 
# Cox model. Because the primary goal was to estimate power under observed trial 
# conditions—not to minimize bias at the cost of further reduced sensitivity—the 
# standard Cox model was used for the main power simulation.

cocaine_clean <- cocaine_data_survival %>%
  select(days_recoded, status_recoded, condition_recoded, Gender) %>%
  filter(
    !is.na(days_recoded),
    !is.na(status_recoded),
    !is.na(condition_recoded),
    !is.na(Gender),
    is.finite(days_recoded),
    is.finite(status_recoded)
  ) %>%
  mutate(
    days_recoded = ifelse(days_recoded == 0, 0.01, days_recoded),
    condition_recoded = factor(condition_recoded, levels = c(0, 1), labels = c("Placebo", "Psilocybin")),
    Gender = factor(as.character(Gender))
  )

firth_model <- coxphf(
  Surv(days_recoded, status_recoded) ~ condition_recoded,
  data = cocaine_clean
)

summary(firth_model)

obs_hr_firth <- exp(coef(firth_model)["condition_recodedPsilocybin"])
event_rate <- mean(cocaine_data_survival$status_recoded)
n <- 40
n_sim <- 1000
alpha <- 0.05

p_vals_firth <- numeric(n_sim)

set.seed(2025)

for (i in 1:n_sim) {
  sim_data <- simulate_cox_data(n = n, hr = obs_hr_firth, event_rate = event_rate)
  sim_data$time <- ifelse(sim_data$time == 0, 0.01, sim_data$time)
  sim_data$treatment <- factor(sim_data$treatment, levels = c(0, 1), labels = c("Placebo", "Psilocybin"))
  
  fit <- tryCatch(
    coxphf(Surv(time, status) ~ treatment, data = sim_data),
    error = function(e) return(NULL)
  )
  
  if (!is.null(fit)) {
    p_vec <- summary(fit)$prob
    if ("treatmentPsilocybin" %in% names(p_vec)) {
      p_vals_firth[i] <- p_vec["treatmentPsilocybin"]
    } else {
      p_vals_firth[i] <- NA
    }
  } else {
    p_vals_firth[i] <- NA
  }
  
  if (i %% 10 == 0) cat("Simulation", i, "complete\n")
}

table(is.na(p_vals_firth))
simulated_power_firth_true <- mean(p_vals_firth < alpha, na.rm = TRUE)
simulated_power_firth_true
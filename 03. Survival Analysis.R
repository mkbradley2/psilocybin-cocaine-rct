# 03. Survival Analysis.r
# Psilocybin for Cocaine Use Disorder RCT - Survival Analysis
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-04-05
#
# Purpose:
#   Conduct survival analysis on time to lapse using Cox models and Kaplan-Meier curves.
#   Assess statistical power post-hoc through simulation.

library(haven)       # Alternate Data Formats
library(tidyverse)   # Tidying
library(survival)    # Survival models
library(survminer)   # KM plots
library(broom)       # Model tidying
library(coxphf)      # Firth’s correction

setwd("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Input")

# Load data

cocaine_data <- read.csv("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output/cocaine_RCT_data.csv")

# ================================================
# Data Cleaning for Modeling
# ================================================

cocaine_data_survival <- cocaine_data %>%
  mutate(across(where(is.labelled), as_factor)) %>%
  mutate(
    Condition = as.character(Condition),
    STATUS_reconstructed = case_when(
      # Has DAYS (time to lapse) = EVENT (STATUS 0)
      !is.na(DAYS) ~ 0,
      # No DAYS (never lapsed) = CENSORED (STATUS 1)
      # This includes both completers AND lost to follow-up
      is.na(DAYS) ~ 1,
      TRUE ~ NA_real_
    ),
    
    # Convert to standard survival coding
    Lapse_Event = case_when(
      STATUS_reconstructed == 0 ~ 1,  # Event (lapsed)
      STATUS_reconstructed == 1 ~ 0,  # Censored (no lapse observed)
      TRUE ~ NA_real_
    ),
    
    condition_recoded = case_when(
      Condition == "Psilocybin" ~ 1,
      Condition == "Placebo" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE) %>%
  mutate(
    days_recoded = case_when(
      STATUS == 1 ~ Integration + 180, # If no lapse, observed through Integration + 180  
      ID_Drug == 1027 ~ Integration + 90 + 1, # No 180 Data
      ID_Drug == 1030 ~ Integration + 90 + 1, # No 180 Data
      ID_Drug == 1037 ~ Integration + 1,      # No F/U Data
      ID_Drug == 1012 ~ 0,                    # No Post-Drug Administration Data
      TRUE ~ DAYS + 1  # Event: use actual time to lapse
    )
  )

# ================================================
# Survival Analysis
# ================================================

km_fit <- survfit(Surv(time = cocaine_data_survival$days_recoded, event = cocaine_data_survival$Lapse_Event) ~ condition_recoded, data = cocaine_data_survival)
logrank_test <- survdiff(Surv(days_recoded, Lapse_Event) ~ condition_recoded, data = cocaine_data_survival)
chisq_stat <- logrank_test$chisq
df <- length(logrank_test$n) - 1
p_val <- 1 - pchisq(chisq_stat, df)

cat("Log-Rank Test:\n")
cat("Chi-square =", round(chisq_stat, 3), "\n")
cat("df =", df, "\n")
cat("p-value =", formatC(p_val, digits = 6, format = "f"), "\n")

summary(km_fit)

# Cox Model to get HR and 95% CI
cox_model <- coxph(Surv(days_recoded, Lapse_Event) ~ condition_recoded, data = cocaine_data_survival)
hr <- summary(cox_model)$coef[1, "exp(coef)"]
hr_ci <- summary(cox_model)$conf.int[1, c("lower .95", "upper .95")]

# Annotated label
hr_label <- paste0("HR = ", round(hr, 2), " (95% CI: ", round(hr_ci[1], 2), "–", round(hr_ci[2], 2), ")\n",
                   "Log-rank p = ", formatC(p_val, digits = 4, format = "f"))

# KM plot with risk table and HR annotation
survival_plot <- ggsurvplot(
  km_fit, 
  data = cocaine_data_survival,
  risk.table = TRUE,
  risk.table.title = "Number at Risk (Censored)",
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.2,
  risk.table.fontsize = 4,
  risk.table.col = "strata",
  pval = FALSE,  # we'll use our own annotation instead
  conf.int = TRUE,
  xlab = "Days",
  ylab = "Proportion Without Lapse (%)",
  surv.scale = "percent",
  legend.title = "",
  legend.labs = c("Placebo", "Psilocybin"),
  legend = "right",
  ggtheme = theme_minimal(base_size = 14),
  palette = c("#F8766D", "#00BFC4")
)

# Add HR and p-value as annotation on plot
# survival_plot$plot <- survival_plot$plot + 
#   annotate("text", x = 5, y = 0.05, hjust = 0, vjust = 0, size = 4.5, label = hr_label)

survival_plot$plot <- survival_plot$plot + 
  annotate(
    "text",
    x = 175,  # Adjust based on max follow-up time
    y = 0.95,  # Close to the top
    hjust = 0,
    size = 4.5,
    label = hr_label
  )

print(survival_plot)

ggsave(
  filename = "cocaine_survival_curve_notitle.png",
  plot = print(survival_plot),
  width = 14,
  height = 9,
  dpi = 300
)

surv_object <- Surv(time = cocaine_data_survival$days_recoded, event = cocaine_data_survival$event)

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

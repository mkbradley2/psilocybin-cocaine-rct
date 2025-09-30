# 04. Fisher's Exact and Logistic Regression.r
# Psilocybin for Cocaine Use Disorder RCT - Fisher's Exact Test, Firth Logistic Regression, and Hedges' g
# Author: Melissa Bradley
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-04-05
#
# Purpose:
#   Conduct Fisher's Exact Test and Firth's logistic regression to assess relapse risk.
#   Estimate effect size using Hedges' g.
#   Perform post-hoc power simulations using Firth's logistic regression.

library(tidyverse)   # Data tidying
library(logistf)     # Firth logistic regression
library(effsize)     # Hedges' g

# Load data
cocaine_data <- read_sav("Peters data testing SPSS.sav")
time_vars <- read_sav("Observation Days.sav")
status_and_days <- read_sav("Otto survival.sav")

# ================================================
# Data Cleaning for Modeling
# ================================================

# Merge survival status and filter to complete cases
cocaine_data_filtered <- cocaine_data %>%
  left_join(select(status_and_days, ID_Sub, STATUS, DAYS), by = "ID_Sub") %>%
  filter(ID_Sub %in% time_vars$ID_Sub, !is.na(Time_Recoded)) %>%
  mutate(
    status_recoded = if_else(STATUS == 0, 1, if_else(STATUS == 1, 0, NA_real_)),
    condition_recoded = if_else(Condition == 2, 1, if_else(Condition == 1, 0, NA_real_))
  ) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE)

# Prepare observation duration variables
specific_columns <- time_vars %>%
  group_by(ID_Sub) %>%
  mutate(days_under_obs = rowSums(across(c(Prescreen, Preparation, Drug, Integration)), na.rm = TRUE)) %>%
  select(ID_Sub, Prescreen, Preparation, Drug, Integration, days_under_obs) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup()

# Merge into final modeling dataset
cocaine_data_log_exact <- cocaine_data_filtered %>%
  left_join(specific_columns, by = "ID_Sub") %>%
  distinct(ID_Sub, .keep_all = TRUE) %>%
  mutate(
    days_recoded = case_when(
      STATUS == 1 ~ Integration,
      ID_Drug %in% c(1027, 1030) ~ Integration + 91,
      ID_Drug == 1037 ~ Integration + 1,
      ID_Drug == 1012 ~ 0,
      TRUE ~ if_else(!is.na(DAYS), DAYS + 1, NA_real_)
    )
  ) %>%
  select(ID_Sub, Gender, Ethnicity, Race, Condition, days_recoded, condition_recoded, status_recoded, STATUS, Integration)

# Fisher's Exact Test
fisher_table <- table(cocaine_data_log_exact$condition_recoded, cocaine_data_log_exact$status_recoded)
fisher.test(fisher_table)

# Firth Logistic Regression - No Covariates
model_no_covariates <- logistf(STATUS ~ condition_recoded, data = cocaine_data_log_exact)
summary(model_no_covariates)

# OR table - No Covariates
or_no_cov <- data.frame(
  Term = names(model_no_covariates$coefficients),
  OR = exp(model_no_covariates$coefficients),
  Lower_CI = exp(model_no_covariates$ci.lower),
  Upper_CI = exp(model_no_covariates$ci.upper),
  p_value = round(model_no_covariates$prob, 8)
) %>%
  mutate(across(2:4, round, 3))

# Firth Logistic Regression - With Covariates
model_with_covariates <- logistf(STATUS ~ condition_recoded + Integration, data = cocaine_data_log_exact)
summary(model_with_covariates)

# OR table - With Covariates (**Covariate Models = Supplemental**)
or_with_cov <- data.frame(
  Term = names(model_with_covariates$coefficients),
  OR = exp(model_with_covariates$coefficients),
  Lower_CI = exp(model_with_covariates$ci.lower),
  Upper_CI = exp(model_with_covariates$ci.upper),
  p_value = round(model_with_covariates$prob, 3)
) %>%
  mutate(across(2:4, round, 3))

print(or_no_cov)
print(or_with_cov)

# Hedges' g - Binary outcome: sustained abstinence
hedges_abstinence_result <- cohen.d(status_recoded ~ factor(condition_recoded), data = cocaine_data_log_exact, pooled = TRUE, hedges.correction = TRUE)
print(hedges_abstinence_result)

# Summary table for Hedges' g (sustained abstinence)
hedges_summary <- tibble(
  Analysis = "Hedges' g (Sustained Abstinence)",
  Term = "Effect Size (Psilocybin vs Placebo)",
  Estimate = hedges_abstinence_result$estimate,
  `Lower CI` = hedges_abstinence_result$conf.int[1],
  `Upper CI` = hedges_abstinence_result$conf.int[2],
  `Std. Error` = sqrt(hedges_abstinence_result$var),
  `Effect Magnitude` = as.character(hedges_abstinence_result$magnitude)
)


print(hedges_summary)

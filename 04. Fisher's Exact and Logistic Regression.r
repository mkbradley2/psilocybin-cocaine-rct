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
library(haven)       # Loading in different data types
library(foreign)     # Loading in different data types p. 2

# Load data
# cocaine_data <- read_sav("Peters data testing SPSS.sav")
# time_vars <- read_sav("Observation Days.sav")
# status_and_days <- read_sav("Otto survival.sav")

cocaine_data <- read.csv("cocaine_RCT_data.csv")

# ================================================
# Data Cleaning for Modeling
# ================================================

cocaine_data_log_exact <- cocaine_data_final %>%
  mutate(
    # Map Time -> ordered integer code
    Time_code = case_when(
      Time == "Prescreening Assessment 1" ~ 1L,
      Time == "Prescreening Assessment 2" ~ 2L,
      Time == "Prescreening Assessment 3" ~ 3L,
      Time == "Prescreening Assessment 4" ~ 4L,
      Time == "Medical Screening"         ~ 5L,
      Time == "Preparation 1"             ~ 6L,
      Time == "Preparation 2"             ~ 7L,
      Time == "Preparation 3"             ~ 8L,
      Time == "Preparation 4"             ~ 9L,
      Time == "Preparation 5"             ~ 10L,
      Time == "MRI 1"                     ~ 11L,
      Time == "Drug Administration"       ~ 12L,
      Time == "Integration/MRI2"          ~ 13L,
      Time == "Follow-up 1"               ~ 14L,
      Time == "Follow-up 2"               ~ 15L,
      Time == "Follow-up 3"               ~ 16L,
      Time == "Follow-up 4"               ~ 17L,
      Time == "90-day Assessment"         ~ 18L,
      Time == "180-day Assessment"        ~ 19L,
      TRUE ~ NA_integer_
    ),
    # Collapse Time_code into broader phases
    Time_Recoded = case_when(
      Time_code == 1L        ~ 1L,   # prescreen 1
      Time_code %in% 2L:6L   ~ 2L,   # prescreen 2..prep 1
      Time_code %in% 7L:12L  ~ 3L,   # prep 2..drug admin
      Time_code %in% 13L:17L ~ 4L,   # integration/followups
      Time_code == 18L       ~ 5L,   # 90-day
      Time_code == 19L       ~ 6L,   # 180-day
      TRUE ~ NA_integer_
    )
  ) %>%
  group_by(ID_Sub) %>%
  # Only fill the fields we actually need
  fill(Time_code, Time_Recoded, STATUS, Condition, ID_Drug, DAYS, Integration,
       .direction = "downup") %>%
  ungroup() %>%
  mutate(
    # status: 1 = event (lapse), 0 = censored
    status_recoded = case_when(
      STATUS == 0 ~ 1L,
      STATUS == 1 ~ 0L,
      TRUE ~ NA_integer_
    ),
    # condition: 1 = Psilocybin, 0 = Placebo (adjust if your coding differs)
    condition_recoded = case_when(
      Condition == "Psilocybin" ~ 1,
      Condition == "Placebo"    ~ 0,
      TRUE ~ NA_real_
    ),
    # recode days; make numeric to avoid type conflicts
    days_recoded = case_when(
      STATUS == 1                    ~ as.numeric(Integration),      # abstinent censored at Integration
      ID_Drug %in% c(1027, 1030)     ~ as.numeric(Integration) + 91, # LTFU at ~90d (+1 for inclusive if desired)
      ID_Drug == 1037                ~ as.numeric(Integration) + 1,
      ID_Drug == 1012                ~ 0,
      !is.na(DAYS)                   ~ as.numeric(DAYS) + 1,
      TRUE                           ~ NA_real_
    )
  ) %>%
  distinct(ID_Sub, .keep_all = TRUE) %>%
  select(
    ID_Sub, Gender, Ethnicity, Race, Condition,
    Prescreen, Preparation, Drug, Integration,
    days_recoded, condition_recoded, status_recoded, STATUS
  )

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

# Firth Logistic Regression - With Integration Days
model_with_covariates <- logistf(STATUS ~ condition_recoded + Integration, data = cocaine_data_log_exact)
summary(model_with_covariates)

# OR table - With Integration Days (**Covariate Models = Supplemental**)
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

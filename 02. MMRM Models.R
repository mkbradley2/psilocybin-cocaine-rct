# 02. MMRM Models.r
# Psilocybin-Cocaine RCT - T-Tests and Mixed Models for Repeated Measures (MMRM)
# Author: Melissa Bradley
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-04-05
#
# Purpose:
#   Fit MMRM models for percent days abstinent from cocaine over time.
#   Conduct t-tests for each timepoint to compare placebo vs. psilocybin.
#   Models include adjusted and unadjusted specifications for pre- and post-
#   treatment periods, with and without condition-by-time interaction terms.
#   Includes estimated marginal means, Hedge's g calculations and 
#   post-hoc power simulations (n=1000) 

library(tidyverse)   # Data Tidying
library(haven)       # Importing data from other software
library(mmrm)        # Mixed Models for Repeated Measures
library(emmeans)     # Estimated marginal means and contrasts (Hedge's g)
library(lme4)        # Linear mixed-effects models for simr power simulation
library(simr)        # Simulation-based power analysis for mixed models

# Load data
cocaine_data <- read_sav("Peters data testing SPSS.sav")
time_vars <- read_sav("Observation Days.sav")
status_and_days <- read_sav("Otto survival.sav")

# ================================================
# Data Cleaning for Modeling
# ================================================

# Specific columns for observation days
specific_columns <- time_vars %>%
  group_by(ID_Sub) %>%
  mutate(days_under_obs = sum(Prescreen, Preparation, Drug, Integration, na.rm = TRUE)) %>%
  select(ID_Sub, days_under_obs) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup()

cocaine_data_t_tests <- cocaine_data %>%
  filter(!is.na(ID_Drug)) %>%
  mutate(
    Cocaine_Variable = case_when(
      Time_Recoded == 1 ~ AB_Per_PreContact,
      Time_Recoded == 2 ~ AB_Per_PreTX,
      Time_Recoded == 3 ~ AB_Per_Preps,
      Time_Recoded == 4 ~ AB_Per_PostDrugTX,
      Time_Recoded %in% c(5, 6) ~ AB_Per_Follow_Ups,
      TRUE ~ NA_real_
    ),
    condition_recoded = case_when(
      Condition == 2 ~ 1,
      Condition == 1 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(
    specific_columns %>%
      group_by(ID_Sub) %>%
      fill(days_under_obs, .direction = "downup") %>%
      summarise(days_under_obs = first(days_under_obs), .groups = "drop"),
    by = "ID_Sub"
  ) %>%
  filter(!is.na(Cocaine_Variable)) %>%
  fill(Gender, .direction = "downup") %>%
  select(ID_Sub, Time_Recoded, Cocaine_Variable, Gender, Condition, condition_recoded, days_under_obs)

valid_post_ids <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(4, 5, 6)) %>%
  pull(ID_Sub) %>%
  unique()

cocaine_data_t_tests <- cocaine_data_t_tests %>%
  filter(ID_Sub %in% valid_post_ids)

cocaine_data_t_tests[sapply(cocaine_data_t_tests, is.infinite)] <- NA

# ================================================
# T-Tests and MMRM Modeling
# ================================================

# T-Tests

print(
  t_test_results <- cocaine_data_t_tests %>%
    group_by(Time_Recoded) %>%
    summarise(
      Mean_Psilocybin = mean(Cocaine_Variable[Condition == 2], na.rm = TRUE),
      Mean_Placebo = mean(Cocaine_Variable[Condition == 1], na.rm = TRUE),
      Mean_Difference = Mean_Psilocybin - Mean_Placebo,
      T_Statistic = t.test(Cocaine_Variable ~ Condition)$statistic,
      P_Value = t.test(Cocaine_Variable ~ Condition)$p.value,
      CI_Lower = t.test(Cocaine_Variable ~ Condition)$conf.int[1],
      CI_Upper = t.test(Cocaine_Variable ~ Condition)$conf.int[2]
    )
)

# MMRM Datasets
cocaine_data_mmrm_full <- cocaine_data_t_tests %>%
  mutate(ID_Sub = as.factor(ID_Sub),
         across(c(condition_recoded, Time_Recoded), as.factor))

cocaine_data_mmrm_pre <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(1, 2, 3)) %>%
  mutate(ID_Sub = as.factor(ID_Sub),
         across(c(condition_recoded, Time_Recoded), as.factor))

cocaine_data_mmrm_post <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(4, 5, 6)) %>%
  mutate(ID_Sub = as.factor(ID_Sub),
         across(c(condition_recoded, Time_Recoded), as.factor))

# Models
model_full_no_inter <- mmrm(Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_full, control = mmrm_control(method = "Satterthwaite"))
model_full_with_inter <- mmrm(Cocaine_Variable ~ condition_recoded * Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_full, control = mmrm_control(method = "Satterthwaite"))
model_pre_no_inter <- mmrm(Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_pre, control = mmrm_control(method = "Satterthwaite"))
model_pre_with_inter <- mmrm(Cocaine_Variable ~ condition_recoded * Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_pre, control = mmrm_control(method = "Satterthwaite"))
model_post_no_inter <- mmrm(Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_post, control = mmrm_control(method = "Satterthwaite"))
model_post_with_inter <- mmrm(Cocaine_Variable ~ condition_recoded * Time_Recoded + days_under_obs + Gender + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_post, control = mmrm_control(method = "Satterthwaite"))
model_post_unadj_no_inter <- mmrm(Cocaine_Variable ~ condition_recoded + Time_Recoded + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_post, control = mmrm_control(method = "Satterthwaite"))
model_post_unadj_with_inter <- mmrm(Cocaine_Variable ~ condition_recoded * Time_Recoded + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_post, control = mmrm_control(method = "Satterthwaite"))
model_pre_unadj_no_inter <- mmrm(Cocaine_Variable ~ condition_recoded + Time_Recoded + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_pre, control = mmrm_control(method = "Satterthwaite"))
model_pre_unadj_with_inter <- mmrm(Cocaine_Variable ~ condition_recoded * Time_Recoded + ar1(Time_Recoded | ID_Sub), data = cocaine_data_mmrm_pre, control = mmrm_control(method = "Satterthwaite"))

extract_mmrm_summary <- function(model, name) {
  summary_df <- as.data.frame(summary(model)$coefficients)
  conf_intervals <- confint(model)
  summary_df <- cbind(summary_df, conf_intervals) %>%
    mutate(
      Model = name,
      Term = rownames(summary_df),
      Estimate = round(Estimate, 3),
      `Std. Error` = round(`Std. Error`, 3),
      `t value` = round(`t value`, 3),
      `Pr(>|t|)` = formatC(`Pr(>|t|)`, digits = 12, format = "f"),
      `CI Lower` = round(`2.5 %`, 3),
      `CI Upper` = round(`97.5 %`, 3)
    ) %>%
    select(Model, Term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`, `CI Lower`, `CI Upper`)
}

mmrm_summaries <- bind_rows(
  extract_mmrm_summary(model_full_no_inter, "Full - No Interaction"),
  extract_mmrm_summary(model_full_with_inter, "Full - With Interaction"),
  extract_mmrm_summary(model_pre_no_inter, "Pre - No Interaction"),
  extract_mmrm_summary(model_pre_with_inter, "Pre - With Interaction"),
  extract_mmrm_summary(model_post_no_inter, "Post - No Interaction"),
  extract_mmrm_summary(model_post_with_inter, "Post - With Interaction"),
  extract_mmrm_summary(model_pre_unadj_no_inter, "Pre Unadjusted - No Interaction"),
  extract_mmrm_summary(model_pre_unadj_with_inter, "Pre Unadjusted - With Interaction"),
  extract_mmrm_summary(model_post_unadj_no_inter, "Post Unadjusted - No Interaction"),
  extract_mmrm_summary(model_post_unadj_with_inter, "Post Unadjusted - With Interaction")
)

print(mmrm_summaries)

emm_results <- emmeans(model_post_with_inter, ~ condition_recoded | Time_Recoded)
contrast_df <- as.data.frame(contrast(emm_results, method = "revpairwise"))
resid_sd <- sqrt(mean(resid(model_post_with_inter)^2, na.rm = TRUE))

hedges_g_results <- contrast_df %>%
  mutate(
    hedges_g = round(estimate / resid_sd, 2),
    estimate = round(estimate, 2),
    SE = round(SE, 2),
    df = round(df, 1),
    t.ratio = round(t.ratio, 2),
    p.value = formatC(p.value, digits = 6, format = "f"),
    resid_sd = round(resid_sd, 2)
  )

print(hedges_g_results)

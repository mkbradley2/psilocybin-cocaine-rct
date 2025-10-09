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
#   Includes estimated marginal means and Hedge's g calculations.

library(tidyverse)   # Data Tidying
library(haven)       # Importing data from other software
library(mmrm)        # Mixed Models for Repeated Measures
library(emmeans)     # Estimated marginal means and contrasts (Hedge's g)
library(lme4)        # Linear mixed-effects models 

# Load data
cocaine_data_final <- read.csv("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output/cocaine_RCT_data.csv")

# ================================================
# Data Cleaning for Modeling
# ================================================

cocaine_data_t_tests <- cocaine_data_final %>%
  filter(!is.na(ID_Drug)) %>%
  mutate(
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
    Time_Recoded = case_when(
      Time_code == 1            ~ 1,
      Time_code %in% 2:6        ~ 2,
      Time_code %in% 7:12       ~ 3,
      Time_code %in% 13:17      ~ 4,
      Time_code == 18           ~ 5,
      Time_code == 19           ~ 6,
      TRUE ~ NA_real_
  ),
    Cocaine_Variable = case_when(
      Time_Recoded == 1        ~ AB_Per_PreContact,
      Time_Recoded == 2        ~ AB_Per_PreTX,
      Time_Recoded == 3        ~ AB_Per_Preps,
      Time_Recoded == 4        ~ AB_Per_PostDrugTX,
      Time_Recoded %in% c(5,6) ~ AB_Per_Follow_Ups,
      TRUE ~ NA_real_
    ),
    condition_recoded = case_when(
      Condition == "Psilocybin" ~ 1,
      Condition == "Placebo"    ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(ID_Sub) %>%
  fill(days_under_obs, .direction = "downup") %>%
  mutate(days_under_obs = first(days_under_obs)) %>%
  ungroup() %>%
  filter(!is.na(Cocaine_Variable)) %>%
  fill(Gender, .direction = "downup") %>%
  fill(AB_Per_PreContact, .direction = "downup") %>%
  select(
    ID_Sub, Time_Recoded, Cocaine_Variable,
    Gender, Condition, condition_recoded, days_under_obs, AB_Per_PreContact
  )


valid_post_ids <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(4, 5, 6)) %>%
  pull(ID_Sub) %>%
  unique()

# cocaine_data_t_tests <- cocaine_data_t_tests %>%
#   filter(ID_Sub %in% valid_post_ids)

# ================================================
# Study Day Variable Creation 
# ================================================

# 1) Lookup stays the same
drug_admin_lookup <- tribble(
  ~ID_Drug_chr, ~DrugAdminDate_str,
  "1001","8/12/15","1002","11/19/15","1003","11/24/15",
  "1005","3/10/16","1006","4/19/16","1007","4/27/16",
  "1008","6/16/16","1009","9/8/16","1010","11/10/16",
  "1011","10/27/16","1012","1/13/17","1013","1/20/17",
  "1014","3/2/17","1015","7/20/17","1016","9/12/17",
  "1017","9/17/18","1004","2/13/19","1018","3/26/19",
  "1023","7/22/19","1019","9/13/19",
  "1020","6/22/20","1021","7/6/20","1022","8/7/20",
  "1024","9/3/20","1025","9/21/20","1026","11/23/20",
  "1027","12/1/20","1028","2/24/21","1029","4/20/21",
  "1030","4/27/21","1031","6/14/21","1032","7/26/21",
  "1033","9/24/21","1034","1/19/22","1035","4/4/22",
  "1036","4/15/22","1037","8/25/22","1038","4/14/23",
  "1039","7/7/23","1040","10/2/23"
) %>%
  mutate(DrugAdminDate = as.Date(DrugAdminDate_str, format = "%m/%d/%y")) %>%
  select(ID_Drug_chr, DrugAdminDate)

# 2) Subject dates (unchanged)
subject_dates <- cocaine_data_final %>%
  transmute(ID_Sub, ID_Drug_chr = as.character(ID_Drug)) %>%
  distinct() %>%
  left_join(drug_admin_lookup, by = "ID_Drug_chr") %>%
  select(ID_Sub, DrugAdminDate)

# 3) Trial start and StudyDay fields (drop old cols first to avoid .x/.y)
trial_start <- min(subject_dates$DrugAdminDate, na.rm = TRUE)

cocaine_data_t_tests <- cocaine_data_t_tests %>%
  select(-any_of(c("DrugAdminDate", "StudyDay", "StudyDay_c"))) %>%  # <-- simple fix
  left_join(subject_dates, by = "ID_Sub") %>%
  mutate(
    StudyDay   = if_else(is.na(DrugAdminDate), NA_integer_,
                         as.integer(DrugAdminDate - trial_start) + 1L),
    StudyDay_c = StudyDay - mean(StudyDay, na.rm = TRUE)
  )

# ================================================
# T-Tests and MMRM Modeling
# ================================================

# T-Tests

print(
  t_test_results <- cocaine_data_t_tests %>%
    group_by(Time_Recoded) %>%
    summarise(
      Mean_Psilocybin = mean(Cocaine_Variable[condition_recoded == 1], na.rm = TRUE),
      Mean_Placebo = mean(Cocaine_Variable[condition_recoded == 0], na.rm = TRUE),
      Mean_Difference = Mean_Psilocybin - Mean_Placebo,
      T_Statistic = t.test(Cocaine_Variable ~ condition_recoded)$statistic,
      P_Value = t.test(Cocaine_Variable ~ condition_recoded)$p.value,
      CI_Lower = t.test(Cocaine_Variable ~ condition_recoded)$conf.int[1],
      CI_Upper = t.test(Cocaine_Variable ~ condition_recoded)$conf.int[2]
    )
)

# MMRM Datasets
cocaine_data_mmrm_full <- cocaine_data_t_tests %>%
  mutate(ID_Sub = as.factor(ID_Sub),
         across(c(condition_recoded, Time_Recoded), as.factor))

cocaine_data_mmrm_post <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(4, 5, 6)) %>%
  mutate(ID_Sub = as.factor(ID_Sub),
         across(c(condition_recoded, Time_Recoded), as.factor))

# -------- FULL --------
model_full_no_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_full,
  control = mmrm_control(method = "Satterthwaite")  # explicit default
)

model_full_with_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded * Time_Recoded + days_under_obs + Gender +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_full,
  control = mmrm_control(method = "Satterthwaite")
)

# -------- POST --------
model_post_no_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)

model_post_with_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded * Time_Recoded + days_under_obs + Gender +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)

# -------- POST (unadjusted) --------
model_post_unadj_no_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded + Time_Recoded +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)

model_post_unadj_with_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded * Time_Recoded +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)
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
  extract_mmrm_summary(model_post_no_inter, "Post - No Interaction"),
  extract_mmrm_summary(model_post_with_inter, "Post - With Interaction"),
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

# Study Day Models

# ================================================
# Adjusted Models (Gender and Study Time)
# ================================================

model_post_adj_no_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded + Time_Recoded + days_under_obs + Gender + StudyDay_c +
    ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)

model_post_adj_with_inter <- mmrm(
  Cocaine_Variable ~ condition_recoded * Time_Recoded + StudyDay_c + StudyDay_c * condition_recoded + AB_Per_PreContact +
  ar1(Time_Recoded | ID_Sub),
  data = cocaine_data_mmrm_post,
  control = mmrm_control(method = "Satterthwaite")
)

mmrm_adj_summaries <- bind_rows(
  extract_mmrm_summary(model_post_adj_no_inter,  "Post Adjusted - No Interaction"),
  extract_mmrm_summary(model_post_adj_with_inter,"Post Adjusted - With Interaction")
)

print(mmrm_adj_summaries)



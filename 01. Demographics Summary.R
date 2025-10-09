# 01. Demographics Summary.r
# Psilocybin for Cocaine Use Disorder RCT - Baseline Demographic and Clinical Summary
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-04-05
#
# Purpose:
#   Generate baseline demographic summaries and between-group comparisons
#   using descriptive statistics and gtsummary tables. 

library(tidyverse)   # Data Tidying
library(haven)       # Importing data from other software
library(gtsummary)   # Summary tables
library(kableExtra)  # Table features (used w/ gtsummary)
library(ggplot2)     # Plotting

# Set Working Directory
setwd("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output")

cocaine_data <- read.csv("cocaine_RCT_data.csv")

# =====================
# Create Codebook
# =====================
codebook <- data.frame(
  Column_Name = names(cocaine_data),
  Description = sapply(1:ncol(cocaine_data), function(i) {
    label <- attr(cocaine_data[[i]], "label")
    if (is.null(label)) NA else label
  }),
  Value_Labels = sapply(1:ncol(cocaine_data), function(i) {
    labels <- attr(cocaine_data[[i]], "labels")
    if (is.null(labels)) NA else paste(names(labels), "=", labels, collapse = "; ")
  })
)

print(codebook)


# =====================
# Check Normality: Q-Q Plots
# =====================
qq_vars <- c("Age", "Age_First_Use", "Age_Reg_Use", "Years_Used", "Quit_Attempts", "Quit_Long")

qq_data <- cocaine_data_joined %>%
  select(Condition, all_of(qq_vars)) %>%
  pivot_longer(cols = -Condition, names_to = "Variable", values_to = "Value") %>%
  filter(!is.na(Value))

qq_plot <- ggplot(qq_data, aes(sample = Value, color = Condition)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Variable, scales = "free", ncol = 2) +
  theme_minimal(base_size = 14) +
  labs(title = "Q-Q Plots of Continuous Variables by Condition") +
  theme(legend.position = "top")

print(qq_plot)

# =====================
# Full Sample Summaries
# =====================
demographic_vars <- cocaine_data %>%
  filter(!is.na(ID_Drug)) %>%
  select(ID_Drug, Gender, Ethnicity, Race, Marital, Age, EduYears, LivSit, Income, Veteran, SexOr, Occup, Degree) %>%
  distinct(ID_Drug, .keep_all = TRUE)

# Categorical Summaries
categorical_summary <- map(demographic_vars, ~table(.x, useNA = "ifany"))
categorical_summary <- categorical_summary[names(categorical_summary) %in% c("Gender", "Ethnicity", "Race", "Marital", "Income", "Veteran", "SexOr", "Occup", "Degree", "LivSit")]

# Continuous Summary
continuous_summary <- demographic_vars %>%
  summarise(
    Age_Mean = mean(Age, na.rm = TRUE),
    Age_SD = sd(Age, na.rm = TRUE),
    Age_Median = median(Age, na.rm = TRUE),
    EduYears_Mean = mean(EduYears, na.rm = TRUE),
    EduYears_SD = sd(EduYears, na.rm = TRUE),
    EduYears_Median = median(EduYears, na.rm = TRUE)
  )

# Output summaries
print("=== Continuous Summary ===")
print(continuous_summary)

print("\n=== Categorical Summaries ===")
walk2(names(categorical_summary), categorical_summary, ~{
  cat("\n=== Summary for", .x, "===\n")
  print(as.data.frame(.y))
})

# =====================
# Table 1: By Condition
# =====================
cocaine_demo <- cocaine_data %>%
  filter(!is.na(ID_Drug), !is.na(Condition)) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE) %>%
  mutate(
    Condition = factor(as.numeric(Condition), levels = c(1, 2), labels = c("Placebo", "Psilocybin")),
    Gender = factor(as.numeric(Gender), levels = c(1, 2), labels = c("Male", "Female")),
    Race_Collapsed = factor(ifelse(as.numeric(Race) == 1, "White", "Black"), levels = c("White", "Black")),
    SexOr = case_when(
      ID_Sub == 89 & ID_Drug == 1023 ~ "Pansexual",
      ID_Sub == 53 & ID_Drug == 1016 ~ "Asexual",
      as.numeric(SexOr) == 1 ~ "Homosexual",
      as.numeric(SexOr) == 2 ~ "Heterosexual",
      as.numeric(SexOr) == 3 ~ "Bisexual",
      as.numeric(SexOr) == 4 ~ "None",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Heterosexual", "Homosexual", "Bisexual", "Asexual", "Pansexual", "None")),
    SexOr_Collapsed = factor(ifelse(SexOr == "Heterosexual", "Heterosexual", "Sexual Minority"), levels = c("Heterosexual", "Sexual Minority")),
    Marital_Status = factor(as.numeric(Marital), levels = c(1, 2, 3, 4), labels = c("Married", "Widowed", "Divorced or Separated", "Single, never married")),
    IncomeGroup = factor(ifelse(as.numeric(Income) <= 3, "< $30,000", "≥ $30,000")),
    PrevPsyUse = factor(ifelse(ID_Drug %in% c(1026, 1038, 1040), "Yes", "No")),
    Treatment = factor(as.numeric(Treatment), levels = c(0, 1), labels = c("No", "Yes")),
    Employment_Status = case_when(
      OccStatus == 1 ~ "Employed",
      OccStatus == 2 ~ "Unemployed",
      OccStatus == 3 ~ "Retired/Disabled",
      OccStatus == 4 ~ "Homemaker",
      OccStatus == 5 ~ "Student",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Employed", "Unemployed", "Retired/Disabled", "Homemaker", "Student")),
    Veteran = factor(as.numeric(Veteran), levels = c(0, 1), labels = c("No", "Yes")),
    LivSit_Collapsed = case_when(
      LivSit %in% c(1, 2) ~ "Housed",
      LivSit == 3 ~ "Therapeutic Housing",
      LivSit == 5 ~ "Homeless",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Housed", "Therapeutic Housing", "Homeless")),
    Edu_Collapsed = case_when(
      Degree == 0 ~ "Less than High School",
      Degree == 1 ~ "High School / GED",
      Degree == 2 ~ "Some College / AA",
      Degree == 3 ~ "Bachelor's Degree",
      Degree %in% c(4, 5) ~ "Graduate or Professional Degree",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Less than High School", "High School / GED", "Some College / AA", "Bachelor's Degree", "Graduate or Professional Degree"))
  )

# Generate Table 1

demographic_table <- cocaine_demo %>%
  select(
    Age, Condition, Gender, Race_Collapsed, SexOr_Collapsed, Marital_Status, IncomeGroup, PrevPsyUse,
    Employment_Status, Veteran, LivSit_Collapsed, Edu_Collapsed,
    Age_First_Use, Age_Reg_Use, Treatment, Quit_Attempts, Quit_Long
  ) %>%
  tbl_summary(
    by = Condition,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous() ~ "{median} ({p25}, {p75})"
    ),
    digits = all_categorical() ~ c(0, 1),
    label = list(
      Gender ~ "Gender",
      Race_Collapsed ~ "Race/Ethnicity",
      SexOr_Collapsed ~ "Sexual Identity",
      Marital_Status ~ "Marital Status",
      IncomeGroup ~ "Household Income",
      PrevPsyUse ~ "Previous Psychedelic Use",
      Employment_Status ~ "Employment Status",
      Veteran ~ "Veteran Status",
      LivSit_Collapsed ~ "Living Situation",
      Edu_Collapsed ~ "Education Level",
      Age_First_Use ~ "Age at First Cocaine Use",
      Age_Reg_Use ~ "Age at Regular Cocaine Use",
      Treatment ~ "Previous Cocaine Treatment",
      Quit_Attempts ~ "Number of Cocaine Quit Attempts",
      Quit_Long ~ "Longest Period of Cocaine Cessation (Days)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_categorical() ~ "fisher.test",
      all_continuous() ~ "wilcox.test"
    )
  ) %>%
  bold_labels()

# Output
print(demographic_table)

# gtsave(as_gt(demographic_table), filename = "results/tables/demographic_table_expanded.rtf")

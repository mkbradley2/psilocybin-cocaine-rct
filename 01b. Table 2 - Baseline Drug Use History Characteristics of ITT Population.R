# 01b. Table 2 - Baseline Drug Use History.r
# Psilocybin-Cocaine RCT - Baseline Drug Use History Characteristics (ITT)
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
#
# Purpose:
#   Replicate Table 2: Baseline Drug Use History Characteristics
#   of the Intention-to-Treat Population

library(tidyverse)
library(gtsummary)

# ============================================================
# 1. Load and prepare data
# ============================================================

cocaine_data <- read.csv("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output/cocaine_RCT_data.csv")

# Time recode map
time_recode_map <- tribble(
  ~Time,                          ~Time_code,
  "Prescreening Assessment 1",    1L,
  "Prescreening Assessment 2",    2L,
  "Prescreening Assessment 3",    3L,
  "Prescreening Assessment 4",    4L,
  "Medical Screening",            5L,
  "Preparation 1",                6L,
  "Preparation 2",                7L,
  "Preparation 3",                8L,
  "Preparation 4",                9L,
  "Preparation 5",                10L,
  "MRI 1",                        11L,
  "Drug Administration",          12L,
  "Integration/MRI2",             13L,
  "Follow-up 1",                  14L,
  "Follow-up 2",                  15L,
  "Follow-up 3",                  16L,
  "Follow-up 4",                  17L,
  "90-day Assessment",            18L,
  "180-day Assessment",           19L
)

base_data <- cocaine_data %>%
  filter(!is.na(ID_Drug)) %>%
  left_join(time_recode_map, by = "Time") %>%
  mutate(
    Time_Recoded = case_when(
      Time_code == 1            ~ 1,
      Time_code %in% 2:6        ~ 2,
      Time_code %in% 7:12       ~ 3,
      Time_code %in% 13:17      ~ 4,
      Time_code == 18           ~ 5,
      Time_code == 19           ~ 6,
      TRUE                      ~ NA_real_
    ),
    Condition = factor(Condition, levels = c("Psilocybin", "Placebo"))
  ) %>%
  filter(!is.na(Time_Recoded))

# ============================================================
# 2. Extract intake (Time_Recoded == 1) data
# ============================================================

intake <- base_data %>%
  group_by(ID_Sub) %>%
  fill(Condition, Gender, Age_First_Use, Age_Reg_Use, Quit_Attempts,
       Quit_Long, Treatment, .direction = "downup") %>%
  ungroup() %>%
  filter(Time_Recoded == 1) %>%
  distinct(ID_Sub, .keep_all = TRUE)

# ============================================================
# 3. Derive Table 2 variables
# ============================================================

table2_data <- intake %>%
  mutate(
    # --- Route of cocaine administration ---
    used_crack   = !is.na(Days_Crack) & Days_Crack > 0,
    used_powder  = !is.na(Days_Cocaine) & Days_Cocaine > 0,
    Route = case_when(
      used_crack & used_powder  ~ "Both smoking and snorting",
      used_crack                ~ "Smoking",
      used_powder               ~ "Snorting",
      TRUE                      ~ NA_character_
    ) %>% factor(levels = c("Smoking", "Snorting", "Both smoking and snorting")),
    
    # --- Cocaine use history (continuous, median/IQR) ---
    Age_First_Use = as.numeric(Age_First_Use),
    Age_Reg_Use   = as.numeric(Age_Reg_Use),
    Quit_Attempts = as.numeric(Quit_Attempts),
    Quit_Long     = as.numeric(Quit_Long),
    
    # --- Lifetime treatment for cocaine (binary) ---
    Treatment = factor(
      case_when(
        Treatment == 1 | Treatment == "Yes" ~ "Yes",
        Treatment == 0 | Treatment == "No"  ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    
    # --- Daily tobacco use (Per_Tob == 100 at intake) ---
    Tobacco_Daily = factor(
      ifelse(!is.na(Per_Tob) & Per_Tob == 100, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    
    # --- Any heavy drinking ---
    Any_Heavy_Drink = factor(
      ifelse(!is.na(Per_Heavy_GenSpec) & Per_Heavy_GenSpec > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    Heavy_Drink_Pct = ifelse(Per_Heavy_GenSpec > 0, Per_Heavy_GenSpec, NA_real_),
    
    # --- Any cannabis use ---
    Any_Cannabis = factor(
      ifelse(!is.na(Per_Can) & Per_Can > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    Cannabis_Pct = ifelse(Per_Can > 0, Per_Can, NA_real_),
    
    # --- Any nonmedical depressant use ---
    Any_Depressant = factor(
      ifelse(!is.na(Per_Dep) & Per_Dep > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    
    # --- Any meth/nonmedical stimulant use ---
    Any_Stimulant = factor(
      ifelse(!is.na(Per_Stim) & Per_Stim > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    
    # --- Any nonmedical opioid use ---
    Any_Opioid = factor(
      ifelse(!is.na(Per_Opi) & Per_Opi > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    
    # --- Any heroin use ---
    Any_Heroin = factor(
      ifelse(!is.na(Per_Her) & Per_Her > 0, 1, 0),
      levels = c(0, 1), labels = c("No", "Yes")
    ),
    
    # --- Lifetime history of hallucinogen use ---
    Prev_Psych_Use = factor(
      ifelse(ID_Drug %in% c(1026, 1038, 1040), "Yes", "No"),
      levels = c("No", "Yes")
    )
  )

# ============================================================
# 4. Verify counts against published Table 2
# ============================================================

cat("\n=== Quick verification ===\n\n")

cat("Route of administration:\n")
table2_data %>% count(Condition, Route) %>% print()

cat("\nAny heavy drinking:\n")
table2_data %>% count(Condition, Any_Heavy_Drink) %>% print()

cat("\nAny cannabis use:\n")
table2_data %>% count(Condition, Any_Cannabis) %>% print()

cat("\nDaily tobacco:\n")
table2_data %>% count(Condition, Tobacco_Daily) %>% print()

cat("\nLifetime hallucinogen use:\n")
table2_data %>% count(Condition, Prev_Psych_Use) %>% print()

# ============================================================
# 5. Generate Table 2 with gtsummary
# ============================================================

table2 <- table2_data %>%
  select(
    Condition,
    Route,
    Age_First_Use, Age_Reg_Use, Quit_Attempts,
    Treatment, Quit_Long,
    Tobacco_Daily,
    Any_Heavy_Drink, Heavy_Drink_Pct,
    Any_Cannabis, Cannabis_Pct,
    Any_Depressant, Any_Stimulant, Any_Opioid, Any_Heroin,
    Prev_Psych_Use
  ) %>%
  tbl_summary(
    by = Condition,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous()  ~ "{median} ({p25}-{p75})"
    ),
    digits = list(
      all_categorical() ~ c(0, 0),
      all_continuous()  ~ 1
    ),
    label = list(
      Route           ~ "Preferred route of cocaine administration",
      Age_First_Use   ~ "Age at first cocaine use, median (IQR), y",
      Age_Reg_Use     ~ "Age of onset of regular cocaine use, median (IQR), y",
      Quit_Attempts   ~ "No. of prior cocaine quit attempts, median (IQR)",
      Treatment       ~ "Lifetime history of treatment for cocaine use",
      Quit_Long       ~ "Longest duration of cocaine abstinence, median (IQR), d",
      Tobacco_Daily   ~ "Daily tobacco use in the past 90 d",
      Any_Heavy_Drink ~ "Any heavy drinking in the past 90 d",
      Heavy_Drink_Pct ~ "Heavy drinking days among those reporting any heavy drinking, median (IQR), %",
      Any_Cannabis    ~ "Any cannabis use in the past 90 d",
      Cannabis_Pct    ~ "Cannabis use days among those reporting any cannabis use, median (IQR), %",
      Any_Depressant  ~ "Any nonmedical depressant use in the past 90 d",
      Any_Stimulant   ~ "Any methamphetamine or nonmedical stimulant use in the past 90 d",
      Any_Opioid      ~ "Any nonmedical opioid use in the past 90 d",
      Any_Heroin      ~ "Any heroin use in the past 90 d",
      Prev_Psych_Use  ~ "Lifetime history of hallucinogen use"
    ),
    missing = "no"
  ) %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Participants, No. (%)**") %>%
  bold_labels()

print(table2)

# ============================================================
# 6. Export
# ============================================================

output_dir <- "C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output/"

# HTML
gtsave(as_gt(table2), file.path(output_dir, "Table2_Baseline_Drug_Use_History.html"))

cat("\nTable 2 exported.\n")
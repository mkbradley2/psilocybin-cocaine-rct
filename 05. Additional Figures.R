# 05. Additional Figures.r
# Psilocybin CUD Trial - Waterfall Plots
# Melissa Bradley
# 01/20/25

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Load Data
status_and_days <- read_sav("Otto survival.sav")
time_vars <- read_sav("Observation Days.sav")
cocaine_data <- read_sav("Peters data testing SPSS.sav")
cocaine_data_master <- read.csv("2024_Data_Filament.csv")

# Analysis Data Prep

time_vars <- time_vars %>%
  group_by(ID_Sub) %>%
  mutate(days_under_obs = sum(Prescreen, Preparation, Drug, Integration, na.rm = TRUE)) %>%
  ungroup()

cocaine_data_joined <- cocaine_data %>%
  left_join(select(status_and_days, ID_Sub, STATUS, DAYS), by = "ID_Sub") %>%
  left_join(select(cocaine_data_master, ID_Sub, AB_Per_DA_to_180), by = "ID_Sub")

cocaine_data_filtered <- cocaine_data_joined %>%
  filter(ID_Sub %in% time_vars$ID_Sub, !is.na(Time_Recoded)) %>%
  mutate(
    status_recoded = case_when(STATUS == 0 ~ 1, STATUS == 1 ~ 0, TRUE ~ NA_real_),
    condition_recoded = case_when(Condition == 2 ~ 1, Condition == 1 ~ 0, TRUE ~ NA_real_)
  ) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE)

specific_columns <- time_vars %>%
  select(ID_Sub, Prescreen, Preparation, Drug, Integration, days_under_obs) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup()

combined_data <- cocaine_data_filtered %>%
  left_join(specific_columns, by = "ID_Sub") %>%
  distinct(ID_Sub, .keep_all = TRUE)

cocaine_data_survival <- combined_data %>%
  mutate(
    days_recoded = case_when(
      STATUS == 1 ~ Integration,
      ID_Drug %in% c(1027, 1030) ~ Integration + 91,
      ID_Drug == 1037 ~ Integration + 1,
      ID_Drug == 1012 ~ 0,
      TRUE ~ if_else(!is.na(DAYS), DAYS + 1, NA_real_)
    )
  ) %>%
  select(ID_Sub, Gender, Ethnicity, Race, Condition, days_recoded, condition_recoded, status_recoded, STATUS, Integration, AB_Per_Preps, AB_Per_PreContact, AB_Per_DA_to_180)

graph_data <- cocaine_data_filtered %>%
  mutate(
    abstinence_baseline = AB_Per_PreContact,
    abstinence_preadmin = AB_Per_Preps,
    abstinence_post = AB_Per_DA_to_180,
    abstinence_change = abstinence_post - abstinence_baseline,
    prep_change = abstinence_post - abstinence_preadmin
  ) %>%
  filter(!is.na(abstinence_baseline) & !is.na(abstinence_post)) %>%
  mutate(Group = if_else(condition_recoded == 1, "Psilocybin", "Placebo"))

# Baseline to Follow-Up Waterfall
waterfall_data_baseline <- graph_data %>%
  arrange(desc(abstinence_change)) %>%
  mutate(Patient = row_number())

waterfall_plot_baseline <- ggplot(waterfall_data_baseline, aes(x = Patient, y = abstinence_change, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Placebo" = "#F8766D", "Psilocybin" = "#00BFC4")) +
  labs(x = "Participants (Ranked by Change in Abstinent Days)", y = "Change in Abstinent Days (%)", fill = "Treatment Group") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5), axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12))

print(waterfall_plot_baseline)
ggsave("waterfall_baseline.png", plot = waterfall_plot_baseline, width = 14, height = 7.5, dpi = 300)

# Prep to Follow-Up Waterfall
waterfall_data_prep <- graph_data %>%
  arrange(desc(prep_change)) %>%
  mutate(Patient = row_number())

waterfall_plot_prep <- ggplot(waterfall_data_prep, aes(x = Patient, y = prep_change, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Placebo" = "#F8766D", "Psilocybin" = "#00BFC4")) +
  labs(x = "Participants (Ranked by Change in Abstinent Days)", y = "Change in Abstinent Days (%)", fill = "Treatment Group") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", size = 14, hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5), axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12))

print(waterfall_plot_prep)
ggsave("waterfall_prep.png", plot = waterfall_plot_prep, width = 14, height = 7.5, dpi = 300)

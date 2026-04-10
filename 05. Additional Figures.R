# 05. Additional Figures.r
# Psilocybin CUD Trial - Waterfall Plots
# Melissa Bradley
# 01/20/25

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Load data
setwd("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Output")
cocaine_data <- read.csv("cocaine_RCT_data.csv")

# ============================================================================
# Data Prep
# ============================================================================

cocaine_data_filtered <- cocaine_data %>%
  filter(!is.na(Time)) %>%
  mutate(
    # STATUS: 0 = Lapsed, 1 = Did not lapse
    # Recode so 1 = event (lapse) for survival modeling
    # STATUS: 0 = Lapsed, 1 = Did not lapse
    # Recode so 1 = event (lapse) for survival modeling
    status_recoded = case_when(
      as.character(STATUS) %in% c("0", "Lapsed")          ~ 1,
      as.character(STATUS) %in% c("1", "Did not lapse")   ~ 0,
      TRUE ~ NA_real_
    ),
    # Condition: 1 = Placebo, 2 = Psilocybin
    # Recode so 1 = Psilocybin, 0 = Placebo
    condition_recoded = case_when(
      as.character(Condition) %in% c("2", "Psilocybin")   ~ 1,
      as.character(Condition) %in% c("1", "Placebo")       ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(ID_Sub) %>%
  fill(everything(), .direction = "downup") %>%
  ungroup() %>%
  distinct(ID_Sub, .keep_all = TRUE)

# ============================================================================
# Waterfall Plot Data
# ============================================================================

graph_data <- cocaine_data_filtered %>%
  mutate(
    abstinence_baseline = AB_Per_PreContact,
    abstinence_preadmin = AB_Per_Preps,
    abstinence_post = AB_Per_DA_to_180,
    abstinence_change = abstinence_post - abstinence_baseline,
    prep_change = abstinence_post - abstinence_preadmin
  ) %>%
  filter(!is.na(abstinence_baseline) & !is.na(abstinence_post)) %>%
  mutate(Group = factor(if_else(condition_recoded == 1, "Psilocybin", "Placebo"),
                        levels = c("Psilocybin", "Placebo")))

# ============================================================================
# Baseline to Follow-Up Waterfall Plot
# ============================================================================

waterfall_data_baseline <- graph_data %>%
  arrange(desc(abstinence_change)) %>%
  mutate(Patient = row_number())

waterfall_plot_baseline <- ggplot(waterfall_data_baseline,
                                  aes(x = Patient, y = abstinence_change, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Psilocybin" = "#00BFC4", "Placebo" = "#F8766D")) +
  labs(
    x = "Participants (Ranked by Change in Abstinent Days)",
    y = "Change in Abstinent Days (%)",
    fill = "Treatment Group"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal() +
  guides(fill = guide_legend(title = NULL)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

print(waterfall_plot_baseline)
ggsave("waterfall_baseline.png", plot = waterfall_plot_baseline, width = 14, height = 7.5, dpi = 300)

# ============================================================================
# Prep to Follow-Up Waterfall Plot
# ============================================================================

waterfall_data_prep <- graph_data %>%
  arrange(desc(prep_change)) %>%
  mutate(Patient = row_number())

waterfall_plot_prep <- ggplot(waterfall_data_prep,
                              aes(x = Patient, y = prep_change, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Psilocybin" = "#00BFC4", "Placebo" = "#F8766D")) +
  labs(
    x = "Participants (Ranked by Change in Abstinent Days)",
    y = "Change in Abstinent Days (%)",
    fill = "Treatment Group"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal() +
  guides(fill = guide_legend(title = NULL)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )

print(waterfall_plot_prep)
ggsave("waterfall_prep.png", plot = waterfall_plot_prep, width = 14, height = 7.5, dpi = 300)

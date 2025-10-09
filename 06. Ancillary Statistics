# 07. Ancillary Statistics
# Psilocybin-Cocaine RCT 
# Lab: Dr. Peter Hendricks, UAB Drug Use and Behavior Lab
# Date: 2025-09-21

library(tidyverse)
library(gt)
library(lubridate)

# Set Working Directory
setwd("C:/Users/mbrad/Desktop/Hendricks Lab/Psilocybin and CUD Clinical Trial (Peter)/Data/Input")

cocaine_data <- read.csv("cocaine_RCT_data.csv")

# Date Since First Participant Drug Administration-----------------

cocaine_data_filled <- cocaine_data_filled %>%
  mutate(Group = recode(Condition, `1` = "Placebo", `2` = "Psilocybin"))

table(is.na(cocaine_data_filled$DrugAdminDate))
table(cocaine_data_filled$Group, useNA = "ifany")
summary(cocaine_data_filled$StudyDay)

with(cocaine_data_filled, t.test(StudyDay ~ Group))

# 1) Group summaries 
studyday_group <- cocaine_data_filled %>%
  filter(!is.na(Group), !is.na(StudyDay)) %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean = mean(StudyDay),
    sd = sd(StudyDay),
    median = median(StudyDay),
    IQR = IQR(StudyDay),
    .groups = "drop"
  ) %>%
  # ensure a consistent order
  arrange(factor(Group, levels = c("Placebo","Psilocybin")))

summary_table <- studyday_group %>%
  mutate(
    across(c(mean, sd, median, IQR), ~round(., 1))
  ) %>%
  gt() %>%
  tab_header(title = "StudyDay by Randomization Group") %>%
  cols_label(
    Group  = "Group",
    n      = "n",
    mean   = "Mean",
    sd     = "SD",
    median = "Median",
    IQR    = "IQR"
  )

# 2) Balance tests

# Q-Q Plots

par(mfrow = c(1,2))  # 2 plots side by side

with(subset(cocaine_data_filled, Group == "Placebo"), {
  qqnorm(StudyDay, main = "Placebo", ylab = "Sample Quantiles")
  qqline(StudyDay, col = "red")
})

with(subset(cocaine_data_filled, Group == "Psilocybin"), {
  qqnorm(StudyDay, main = "Psilocybin", ylab = "Sample Quantiles")
  qqline(StudyDay, col = "blue")
})

par(mfrow = c(1,1))  # reset layout

ttest  <- t.test(StudyDay ~ Group, data = cocaine_data_filled)
wilcox <- wilcox.test(StudyDay ~ Group, data = cocaine_data_filled)

# mean diff = Psilocybin - Placebo
mean_diff <- unname(diff(ttest$estimate))  # (Psilo - Placebo)

# standardized mean difference (using pooled SD)
sd_placebo    <- studyday_group$sd[studyday_group$Group == "Placebo"]
sd_psilocybin <- studyday_group$sd[studyday_group$Group == "Psilocybin"]
sd_pooled <- sqrt((sd_placebo^2 + sd_psilocybin^2) / 2)
smd <- mean_diff / sd_pooled

tests_table <- tibble::tibble(
  Test      = c("Welch t-test", "Wilcoxon rank-sum"),
  `Mean diff (Psilo - Plac)` = c(round(mean_diff, 1), NA),
  Statistic = c(round(as.numeric(ttest$statistic), 2), round(as.numeric(wilcox$statistic), 0)),
  `p-value` = c(signif(ttest$p.value, 3), signif(wilcox$p.value, 3))
) %>%
  gt() %>%
  tab_header(title = "Balance Tests: StudyDay by Group") %>%
  tab_source_note(md(paste0("Standardized mean difference (SMD): **",
                            round(smd, 3), "**")))

# View/export:
summary_table
tests_table

# gtsave(summary_table, "studyday_group_summary.html")
# gtsave(tests_table,  "studyday_balance_tests.html")

# Days Under Observation----------------

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
  select(ID_Sub, Time_Recoded, Cocaine_Variable, Gender, Condition, condition_recoded, days_under_obs, AB_Per_PreContact)

valid_post_ids <- cocaine_data_t_tests %>%
  filter(Time_Recoded %in% c(4, 5, 6)) %>%
  pull(ID_Sub) %>%
  unique()

cocaine_data_t_tests <- cocaine_data_t_tests %>%
  filter(ID_Sub %in% valid_post_ids)

cocaine_data_t_tests[sapply(cocaine_data_t_tests, is.infinite)] <- NA

# Build one record per participant with raw period totals
days_by_id <- time_vars %>%
  group_by(ID_Sub) %>%
  summarise(
    Days_Prescreen   = sum(Prescreen,   na.rm = TRUE),
    Days_Preparation = sum(Preparation, na.rm = TRUE),
    Days_Integration = sum(Integration, na.rm = TRUE),
    .groups = "drop"
  )

# Add group labels (Condition: 1 = Placebo, 2 = Psilocybin)
group_map <- cocaine_data %>%
  select(ID_Sub, Condition) %>%
  distinct() %>%
  mutate(Group = factor(Condition, levels = c(1, 2),
                        labels = c("Placebo", "Psilocybin")))

days_grouped <- days_by_id %>%
  inner_join(group_map, by = "ID_Sub")

# Quick helper to compute t-test + Hedges’ g for a given variable
summ_period <- function(var) {
  t_res <- t.test(reformulate("Group", response = var), data = days_grouped)
  g_res <- effsize::cohen.d(reformulate("Group", response = var),
                            data = days_grouped, hedges.correction = TRUE)
  desc  <- days_grouped %>%
    group_by(Group) %>%
    summarise(
      mean = mean(.data[[var]], na.rm = TRUE),
      sd   = sd(.data[[var]], na.rm = TRUE),
      n    = n(),
      .groups = "drop"
    )
  list(t = t_res, g = g_res, desc = desc)
}

# Run for each period
res_prescreen   <- summ_period("Days_Prescreen")
res_preparation <- summ_period("Days_Preparation")
res_integration <- summ_period("Days_Integration")

# View descriptive stats and tests
res_prescreen$desc;   res_prescreen$t;   res_prescreen$g
res_preparation$desc; res_preparation$t; res_preparation$g
res_integration$desc; res_integration$t; res_integration$g

# CEQ Credibility and Expectancy----

#--- 1) Average the three items per subscale per participant per timepoint ---
ceq_avg <- cocaine_data %>%
  select(ID_Sub, Condition, Time,
         CQ_Cred1, CQ_Cred2, CQ_Cred3,
         CQ_Exp1,  CQ_Exp2,  CQ_Exp3) %>%
  pivot_longer(
    cols = starts_with("CQ_"),
    names_to = c("Scale", "Item"),
    names_pattern = "CQ_(Cred|Exp)([123])",
    values_to = "score"
  ) %>%
  group_by(ID_Sub, Condition, Time, Scale) %>%
  summarise(score_mean = mean(score, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(score_mean))

#--- 2) Q-Q plots to assess normality BEFORE choosing test ---
qq_data <- ceq_avg %>%
  filter(!is.na(Condition)) %>%
  mutate(
    Condition_lbl = case_when(
      Condition %in% c(1, "1", "Psilocybin") ~ "Psilocybin",
      Condition %in% c(0, "0", "Placebo")    ~ "Placebo",
      TRUE ~ as.character(Condition)
    )
  )

qq_plot <- ggplot(qq_data, aes(sample = score_mean, color = Condition_lbl)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(Scale ~ Time, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(title = "Q-Q Plots of CEQ Scale Scores by Timepoint and Condition") +
  theme(legend.position = "top")

print(qq_plot)

#--- 3) Effect size helper: rank-biserial correlation ---
rank_biserial <- function(x, y) {
  res <- wilcox.test(x, y, exact = FALSE, correct = TRUE)
  rbc <- (res$statistic / (length(x) * length(y))) - 0.5
  as.numeric(rbc)
}

#--- 4) Run Wilcoxon tests per Scale × Time ---
ceq_results_wilcox <- ceq_avg %>%
  filter(!is.na(Condition)) %>%
  group_by(Scale, Time) %>%
  group_modify(function(df, key) {
    if (n_distinct(df$Condition) < 2) return(tibble())
    
    x <- df %>% filter(Condition %in% c(1, "Psilocybin")) %>% pull(score_mean)
    y <- df %>% filter(Condition %in% c(0, "Placebo"))    %>% pull(score_mean)
    
    summ <- df %>%
      group_by(Condition) %>%
      summarise(
        n  = n(),
        M  = mean(score_mean, na.rm = TRUE),
        SD = sd(score_mean, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(names_from = Condition, values_from = c(n, M, SD))
    
    ww <- wilcox.test(x, y, exact = FALSE, correct = TRUE)
    
    tibble(
      n_Psilocybin = summ$n_1 %||% summ$n_Psilocybin,
      n_Placebo    = summ$n_0 %||% summ$n_Placebo,
      M_Psilocybin = summ$M_1 %||% summ$M_Psilocybin,
      SD_Psilocybin= summ$SD_1 %||% summ$SD_Psilocybin,
      M_Placebo    = summ$M_0 %||% summ$M_Placebo,
      SD_Placebo   = summ$SD_0 %||% summ$SD_Placebo,
      W            = unname(ww$statistic),
      p            = formatC(ww$p.value, format = "f", digits = 4),
      RankBiserial = round(rank_biserial(x, y), 2)
    )
  }) %>%
  ungroup() %>%
  mutate(
    across(c(M_Psilocybin, SD_Psilocybin, M_Placebo, SD_Placebo), ~ round(., 2))
  ) %>%
  arrange(Scale, Time)

#--- 5) View test results ---
print(ceq_results_wilcox)



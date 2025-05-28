# Psilocybin-Assisted Therapy for Cocaine Use Disorder

This repository contains statistical analysis code for a randomized controlled trial comparing psilocybin-assisted therapy to placebo in individuals with cocaine use disorder.

## Overview

- **Primary Outcomes:**  
  - Percentage of days abstinent from cocaine  
  - Time to first cocaine lapse  
  - Sustained abstinence from cocaine  

- **Statistical Models and Tests:**  
  - T-Tests
  - MMRM (Mixed Models for Repeated Measures)  
  - Kaplan-Meier curve and Cox proportional hazards (Survival Analysis)
  - Fisher’s Exact Test 
  - Firth’s penalized logistic regression   
  - Hedges’ g effect sizes  
  - Post-hoc power simulations for all models  

## Scripts

- **01. Demographics Summary.R**  
  Baseline demographic and clinical characteristics; between-group comparisons using `gtsummary`.

- **02. MMRM Models.R**  
  MMRM analysis of <i>Percentage of days abstinent from cocaine</i> (Primary Outcome 1), including adjusted/unadjusted models, time-by-condition interactions, and Hedge’s g effect size estimation. Includes post-hoc power simulation using `simr`.

- **03. Survival Analysis.R**  
  Kaplan-Meier curves and log-rank tests, Cox models for <i> time to first cocaine lapse</i> (Primary Outcome 2), and survival-based post-hoc power simulations. Includes optional models using Firth's correction.

- **04. Fisher's Exact and Logistic Regression.R**  
  Fisher's Exact Test, Firth logistic regression for <i>Sustained abstinence from cocaine</i> (Primary Outcome 3), Hedges’ g calculation for group effect size, and simulation-based power estimation for binary outcomes.

- **05. Additional Figures.R**
  Waterfall plots for each participant's change in abstinent days (%) by condition (Placebo or Psilocybin) 1) Baseline to Follow-Up and 2) Prepartion to Follow-Up.

- **06. [cud_trial_data_cleaning_notes.md](cud_trial_data_cleaning_notes.md)**  
  Details on how study timepoints were mapped to analysis periods and how survival and abstinence data were prepared, including known issues and resolutions related to missing `DAYS` variables.

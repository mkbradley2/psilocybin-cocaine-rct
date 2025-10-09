# Psilocybin-assisted Psychotherapy in the Treatment of Cocaine Use Disorder: A Randomized Controlled Trial
### Created by: Melissa Bradley
### Version: 05/27/2025

## Background
This repository contains statistical analysis code for a randomized controlled trial of psilocybin-assisted psychotherapy in the treatment of cocaine use disorder.

## Overview

- **Primary Outcomes:**  
  - Percentage of days abstinent from cocaine  
  - Complete abstinence from cocaine
  - Time to first cocaine lapse after drug administration  

- **Statistical Models and Tests:**  
  - Wilcoxon rank-sum tests and independent samples t-tests
  - Mixed Models for Repeated Measures (MMRMs)  
  - Kaplan-Meier curve and Cox proportional hazards (Survival Analysis)
  - Fisher’s Exact Test 
  - Firth’s penalized logistic regression   
  - Hedges’ g effect sizes  
  

## Scripts

- **01. Demographics Summary.R**  
  Baseline demographic and clinical characteristics; between-group comparisons using `gtsummary`.

- **02. MMRM Models.R**  
  MMRM analysis of <i>Percentage of days abstinent from cocaine</i> (Primary Outcome 1), including adjusted/unadjusted models, time-by-condition interactions, and Hedge’s g effect size estimation. 

- **03. Survival Analysis.R**  
  Kaplan-Meier curves and log-rank tests, Cox models for <i> time to first cocaine lapse</i> (Primary Outcome 2), and survival-based post-hoc power simulations. Includes optional models using Firth's correction.

- **04. Fisher's Exact and Logistic Regression.R**  
  Fisher's Exact Test, Firth logistic regression for <i>Sustained abstinence from cocaine</i> (Primary Outcome 3), Hedges’ g calculation for group effect size, and simulation-based power estimation for binary outcomes.

- **05. Additional Figures.R** <br>
  Waterfall plots for each participant's change in abstinent days (%) by condition (Placebo or Psilocybin) 1) Baseline to Follow-Up and 2) Prepartion to Follow-Up.

- **06. Ancillary Statistics.R** <br>
  Additional statistical tests for manuscript.

- **07. [cud_trial_data_cleaning_notes.md](cud_trial_data_cleaning_notes.md)**  
  Details on how study timepoints were mapped to analysis periods and how survival and abstinence data were prepared, including known issues and resolutions related to missing `DAYS` variables.

- **08. [psilocybin_cocaine_RCT_codebook.xlsx](psilocybin_cocaine_RCT_codebook.xlsx)** <br>
  Codebook with variable names, descriptions and values (dataset available upon request from PI on Synapse).

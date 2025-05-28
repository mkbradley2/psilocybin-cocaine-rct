# Psilocybin-Assisted Therapy for Cocaine Use Disorder

This repository contains statistical analysis code for a randomized controlled trial comparing psilocybin-assisted therapy to placebo in individuals with cocaine use disorder.

## Overview

- **Primary Outcomes:**  
  - Percentage of days abstinent from cocaine  
  - Time to first cocaine lapse  
  - Sustained abstinence from cocaine  

- **Statistical Models and Tests:**  
  - MMRM (Mixed Models for Repeated Measures)  
  - Cox proportional hazards survival analysis  
  - Firth’s penalized logistic regression  
  - Fisher’s Exact Test  
  - Hedges’ g effect size  
  - Post-hoc power simulations  

## Scripts

- **01_demographics_summary.R**  
  Baseline demographic and clinical characteristics; between-group comparisons using `gtsummary`.

- **02_mmrm_models.R**  
  MMRM analysis of <i>Percentage of days abstinent from cocaine</i> (Primary Outcome 1), including adjusted/unadjusted models, time-by-condition interactions, and Hedge’s g effect size estimation. Includes post-hoc power simulation using `simr`.

- **03_survival_analysis.R**  
  Kaplan-Meier curves and log-rank tests, Cox models for <i> time to first cocaine lapse</i> (Primary Outcome 2), and survival-based post-hoc power simulations. Includes optional models using Firth's correction.

- **04_logistic_firth_analysis.R**  
  Fisher's Exact Test, Firth logistic regression for <i>Sustained abstinence from cocaine</i> (Primary Outcome 3), Hedges’ g calculation for group effect size, and simulation-based power estimation for binary outcomes.


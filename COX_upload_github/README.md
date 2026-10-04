# Time-Dependent Cox Analyses of Vaccination and Incident Dementia / Alzheimer's Disease

## Overview

This repository contains the R code used for the time-dependent Cox proportional
hazards analyses of vaccination and incident dementia / Alzheimer's disease (AD)
in a prospective cohort.

Three vaccines are evaluated as time-dependent exposures:

- Influenza vaccine (流感)
- Herpes zoster vaccine (带疱)
- Pneumococcal vaccine (肺炎)

Key design features:

- **14-day lag** after vaccination before exposure is considered (`lag_days = 14`)
- **365-day protection window** for the influenza vaccine (`protection_days = 365`)
- Time-dependent vaccine status built with `survival::tmerge`
- Time-dependent disease covariates (hypertension, coronary heart disease,
  stroke, diabetes, thyroid disease, COPD, chronic kidney disease, Parkinson's
  disease, depression) built from baseline prevalence plus follow-up onset
  records, with last-observation-carried-forward (LOCF) imputation
- For the vaccine-plus-disease models, the results are reported as **three Cox
  models**: a base model and two sensitivity analyses (see Section 5 of the code)

## Repository structure

```
COX_upload_github/
├── README.md            # this file
└── COX_upload_github.R  # main analysis script
```

## Requirements

R (>= 4.x) with the following packages:

```r
install.packages(c("dplyr", "ggplot2", "tidyr", "tidyverse", "lubridate",
                   "survival", "purrr", "broom", "tibble", "cmprsk", "data.table"))
```

## Required upstream data objects

The script expects the following objects to already exist in the R session
before it is sourced. The individual-level data are **not** provided, due to
data management regulations.

| Object                | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `final_data_merged_2` | Vaccination records: `高峰号` (subject ID), `接种时间` (vaccination date), `疫苗` (vaccine type: 流感/带疱/肺炎), `days_from_enrollment` |
| `data_59701`          | Analysis cohort: baseline covariates, censoring dates (`dementia_censor_date`, `AD_censor_date`, `death_time`), and outcomes (`dementia_event`, `AD_event`) |
| `plot_data_dedup`     | Long-format disease onset records: `高峰号`, `days_from_enrollment`, `disease` |
| `data0_67572`         | Additional variables, e.g., continuous age (`年龄`)                          |

## Analyses included

1. **Influenza + herpes zoster** — dementia and AD outcomes (Section 3)
2. **Influenza + herpes zoster + pneumococcal** — dementia and AD outcomes (Section 4)
3. **Vaccines + time-dependent diseases** — dementia and AD outcomes (Section 5),
   each reported as:
   - Base model (primary analysis)
   - Sensitivity analysis 1: additionally adjusted for coronary heart disease
     and stroke (`data_CHD`, `data_stroke`)
   - Sensitivity analysis 2: additionally adjusted for Parkinson's disease,
     depression, and family history of PD / AD / depression
4. **Stratified analyses** (Section 6): by four age groups (40-49, 50-59, 60-69,
   70-79); by age <60 / ≥60; by anxiety/depression status among ≥60
5. **Competing-risk analysis for death** (Section 7)
6. **Interaction models** (Section 8): pneumococcal vaccination × age group,
   pneumococcal × influenza, and vaccine × BMI with a unified output table

## How to run

1. Load the required data objects into an R session (see above).
2. Run the script either line by line or with:

   ```r
   source("COX_upload_github.R")
   ```

3. Model summaries are printed to the console. Stratified and interaction
   results are stored in `final_results`, `final_resultsde60`, and
   `result_table`.

## Notes

- Chinese characters in the code are **data column names and category values**
  (e.g., `高峰号`, `是否患高血压`, `"流感"`, `"≥60岁"`). They are kept as-is
  because they must match the original data.
- The script is shared as supplementary material for academic use. Please cite
  the accompanying article when using this code.

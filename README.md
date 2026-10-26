# Immuno-metabolic depression approximation in NHANES

R code, data and results for "Prevalence and Metabolic Correlates of a Questionnaire-Based Approximation of Immuno-Metabolic Depression in United States Adults" (NHANES 2015-2016, 2017-2018 and 2021-2023).

## Folders

- `R/` the analysis scripts, numbered in the order they run.
- `data/raw/` the NHANES files used (public CDC data, saved as .rds).
- `output/tables/` the tables as CSV files, and `output/figures/` the figures as PNG and PDF.

## Scripts

| Script | What it does |
|---|---|
| `R/00_packages.R` | Loads the packages (survey, tidyverse, mice, gtsummary, ggplot2 and others) |
| `R/00_config.R` | Thresholds (PHQ-9 >= 10, four-item symptom score >= 6, CRP > 1 mg/L), CRP bridging constants, covariate list, colours |
| `R/01_download.R` | Downloads the NHANES files; skips every file already in `data/raw/` |
| `R/02_clean.R` | Cleans and merges the cycles, defines depression, the symptom score, CRP and the covariates, builds the survey design |
| `R/02b_dag.R` | Draws the directed acyclic graph (Figure S1) |
| `R/03_analysis.R` | Prevalence, the criterion overlap test, metabolic comparisons, sociodemographic models and diagnostics |
| `R/05_sensitivity.R` | Sensitivity analyses (CRP cut-offs, symptom definitions, strata, cycles) |
| `R/06_sensitivity_ipcw_mi.R` | Inclusion weighting and multiple imputation |
| `R/04_tables_figures.R` | Tables 1 and 2, the supplementary tables, Figures 1 and 2 and the flowchart |

## How to run

Open R in this folder and run the scripts in the order above, for example `source("R/02_clean.R")`. Each script reads from `data/` and writes to `data/processed/` and `output/`. A full run takes about two minutes. Package versions are listed in the supplementary material of the paper.

## Licence

Code: MIT. NHANES data: public domain (CDC/NCHS).

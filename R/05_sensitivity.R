source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in   <- here::here("data", "processed", "analytic_dataset.rds")
.self <- here::here("R", "05_sensitivity.R")
.out  <- here::here("data", "processed", "sensitivity_results.rds")
if (file.exists(.out) &&
    file.mtime(.out) > file.mtime(.in) &&
    file.mtime(.out) > file.mtime(.self)) {
  message("  Up to date, skipping")
} else {

df <- readRDS(here::here("data", "processed", "analytic_dataset.rds"))
full_df <- readRDS(here::here("data", "processed", "full_dataset.rds"))

stopifnot(all(model_covars %in% names(df)))

run_imd_analysis <- function(data, label, design = NULL, strat_var = NULL,
                             covars = model_covars) {
  message("  Sensitivity: ", label)
  # Drop the stratifier from covariates and droplevels (single-level factors error)
  covars_use <- setdiff(covars, strat_var)
  if (!is.null(strat_var)) data <- droplevels(data)

  des <- if (!is.null(design)) design
         else svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                        weights = ~wt_mec, nest = TRUE, data = data)
  des_dep <- subset(des, depression == 1)

  n_dep <- sum(data$depression == 1, na.rm = TRUE)
  n_imd <- sum(data$imd_binary == 1, na.rm = TRUE)

  prev <- tryCatch(svymean(~imd_binary, des_dep, na.rm = TRUE),
                   error = function(e) NULL)
  pred <- tryCatch({
    f <- paste("imd_binary ~", paste(covars_use, collapse = " + "))
    svyglm(as.formula(f), design = des_dep, family = quasibinomial())
  }, error = function(e) NULL)

  list(
    label   = label,
    n_total = nrow(data),
    n_dep   = n_dep,
    n_imd   = n_imd,
    prevalence = prev,
    design_df  = survey::degf(des_dep),
    predictors = pred
  )
}

sensitivity <- list()

# ============================================================
# S1: Strict CRP > 3 mg/L
# ============================================================
df_s1 <- df |> mutate(imd_binary = as.integer(dep_subtype_strict == "IMD"))
sensitivity$s1_strict_crp <- run_imd_analysis(df_s1, "CRP > 3 mg/L (strict)")

# ============================================================
# S2: Include acute CRP > 10 (no exclusion)
# ============================================================
# Re-derive from full dataset using crp_raw (preserves acute values)
df_s2 <- full_df |>
  filter(!is.na(phq9_total), !is.na(RIDAGEYR), RIDAGEYR >= 20,
         !is.na(crp_raw),
         if_all(all_of(c("sex", "race_ethnicity", "education", "INDFMPIR")),
                ~ !is.na(.x))) |>
  mutate(
    dep_subtype_s2 = factor(case_when(
      phq9_total < depression_cutoff ~ "No depression",
      aes_score >= aes_threshold & crp_raw > crp_cutoff_primary ~ "IMD",
      TRUE ~ "Non-IMD depression"
    ), c("No depression", "Non-IMD depression", "IMD")),
    imd_binary = as.integer(dep_subtype_s2 == "IMD")
  )
sensitivity$s2_include_acute <- run_imd_analysis(df_s2, "Include acute CRP > 10")

# ============================================================
# S3: Age-stratified (subset full design for correct variance)
# ============================================================
full_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = ~wt_mec, nest = TRUE, data = df)

for (ag in levels(df$age_group)) {
  design_ag <- subset(full_design, age_group == ag)
  df_ag <- df |> filter(age_group == ag)
  key <- paste0("s3_age_", gsub("[^a-z0-9]", "", tolower(ag)))
  sensitivity[[key]] <- run_imd_analysis(df_ag, paste("Age:", ag),
                                         design = design_ag)
}

# ============================================================
# S4: Sex-stratified (subset full design for correct variance)
# ============================================================
for (sx in levels(df$sex)) {
  design_sx <- subset(full_design, sex == sx)
  df_sx <- df |> filter(sex == sx)
  key <- paste0("s4_sex_", tolower(sx))
  sensitivity[[key]] <- run_imd_analysis(df_sx, paste("Sex:", sx),
                                         design = design_sx,
                                         strat_var = "sex")
}

# ============================================================
# S5: Exclude 2021-2023 (test bridging equation impact)
# ============================================================
df_s5 <- df |>
  filter(cycle != "2021-2023") |>
  mutate(wt_mec = wt_blood_orig / 2)  # 2 cycles remaining
sensitivity$s5_no_2021 <- run_imd_analysis(df_s5, "Exclude 2021-2023")

# ============================================================
# S6: Restricted AES (DPQ040 + DPQ080 only)
# ============================================================
aes_restricted <- c("DPQ040", "DPQ080")
if (all(aes_restricted %in% names(df))) {
  df_s6 <- df |>
    mutate(
      aes_restricted = rowSums(pick(all_of(aes_restricted)), na.rm = FALSE),
      dep_subtype_s6 = factor(case_when(
        phq9_total < depression_cutoff ~ "No depression",
        aes_restricted >= 3 & crp > crp_cutoff_primary ~ "IMD",
        TRUE ~ "Non-IMD depression"
      ), c("No depression", "Non-IMD depression", "IMD")),
      imd_binary = as.integer(dep_subtype_s6 == "IMD")
    )
  sensitivity$s6_restricted_aes <- run_imd_analysis(
    df_s6, "Restricted AES (fatigue + psychomotor change, threshold >= 3)")
}

# ============================================================
# S7: AES threshold ≥ 7 (one above Lamers 2024 primary)
# ============================================================
# Primary cut-off stays at >= 6 (Lamers 2024); 7 and 8 are sensitivity analyses.
df_s7 <- df |>
  mutate(
    dep_subtype_s7 = factor(case_when(
      phq9_total < depression_cutoff ~ "No depression",
      aes_score >= 7 & crp > crp_cutoff_primary ~ "IMD",
      TRUE ~ "Non-IMD depression"
    ), c("No depression", "Non-IMD depression", "IMD")),
    imd_binary = as.integer(dep_subtype_s7 == "IMD")
  )
sensitivity$s7_aes_ge7 <- run_imd_analysis(df_s7, "AES threshold >= 7")

# ============================================================
# S8: AES threshold ≥ 8 (two above Lamers 2024 primary)
# ============================================================
df_s8 <- df |>
  mutate(
    dep_subtype_s8 = factor(case_when(
      phq9_total < depression_cutoff ~ "No depression",
      aes_score >= 8 & crp > crp_cutoff_primary ~ "IMD",
      TRUE ~ "Non-IMD depression"
    ), c("No depression", "Non-IMD depression", "IMD")),
    imd_binary = as.integer(dep_subtype_s8 == "IMD")
  )
sensitivity$s8_aes_ge8 <- run_imd_analysis(df_s8, "AES threshold >= 8")

# ============================================================
# S9-S10: Fatigue-only symptom criterion (the one PHQ-9 AES item whose
# direction is unambiguous): DPQ040 >= 2 and DPQ040 = 3, each with CRP > 1
# ============================================================
fatigue_def <- function(data, cut) {
  data |>
    mutate(
      dep_subtype_f = factor(case_when(
        phq9_total < depression_cutoff ~ "No depression",
        DPQ040 >= cut & crp > crp_cutoff_primary ~ "IMD",
        TRUE ~ "Non-IMD depression"
      ), c("No depression", "Non-IMD depression", "IMD")),
      imd_binary = as.integer(dep_subtype_f == "IMD")
    )
}
sensitivity$s9_fatigue_ge2 <- run_imd_analysis(
  fatigue_def(df, 2), "Fatigue item >= 2 (more than half the days) + CRP > 1")
sensitivity$s10_fatigue_3 <- run_imd_analysis(
  fatigue_def(df, 3), "Fatigue item = 3 (nearly every day) + CRP > 1")

# ============================================================
# S12: 2015-2016 CRP on its native DxC 660i scale (no NCHS bridging), as in
# the first submission; 2021-2023 bridging retained
# ============================================================
df_s12 <- full_df |>
  mutate(crp_s12 = if_else(cycle == "2015-2016", crp_native, crp_raw),
         crp_s12 = if_else(crp_s12 > crp_acute_cutoff, NA_real_, crp_s12)) |>
  filter(!is.na(phq9_total), !is.na(crp_s12),
         if_all(all_of(c("sex", "race_ethnicity", "education", "INDFMPIR")),
                ~ !is.na(.x))) |>
  mutate(
    dep_subtype_s12 = factor(case_when(
      phq9_total < depression_cutoff ~ "No depression",
      aes_score >= aes_threshold & crp_s12 > crp_cutoff_primary ~ "IMD",
      TRUE ~ "Non-IMD depression"
    ), c("No depression", "Non-IMD depression", "IMD")),
    imd_binary = as.integer(dep_subtype_s12 == "IMD")
  )
sensitivity$s12_unbridged_2015 <- run_imd_analysis(
  df_s12, "2015-2016 CRP not bridged to the Cobas 6000 scale")

# ============================================================
# S13: three-item symptom score without the psychomotor item
# (sleep + fatigue + appetite, range 0-9, cut-off 5 = half the range, as
# 6 of 12 for the four-item score)
# ============================================================
df_s13 <- df |>
  mutate(
    sym3 = DPQ030 + DPQ040 + DPQ050,
    dep_subtype_s13 = factor(case_when(
      phq9_total < depression_cutoff ~ "No depression",
      sym3 >= 5 & crp > crp_cutoff_primary ~ "IMD",
      TRUE ~ "Non-IMD depression"
    ), c("No depression", "Non-IMD depression", "IMD")),
    imd_binary = as.integer(dep_subtype_s13 == "IMD")
  )
sensitivity$s13_three_item <- run_imd_analysis(
  df_s13, "Three symptom items without psychomotor change, cut-off 5")

# ============================================================
# S11: Model 3 additionally adjusted for leisure-time physical activity
# (harmonised across cycles; instrument discontinuity in 2021-2023)
# ============================================================
sensitivity$s11_pa_adjusted <- run_imd_analysis(
  df, "Additionally adjusted for leisure-time physical activity",
  covars = model_covars_pa)

# --- Summary table ---
sens_summary <- map_dfr(sensitivity, \(s) {
  prev_est <- if (!is.null(s$prevalence)) round(100 * coef(s$prevalence), 1) else NA
  prev_ci <- if (!is.null(s$prevalence)) {
    ci <- confint(s$prevalence, df = s$design_df)
    paste0("(", round(100 * ci[1], 1), "-", round(100 * ci[2], 1), ")")
  } else NA
  tibble(
    analysis   = s$label,
    n_total    = s$n_total,
    n_depressed = s$n_dep,
    n_imd      = s$n_imd,
    imd_prev   = prev_est,
    imd_ci     = prev_ci
  )
})

write_csv(sens_summary, here::here("output", "tables", "sensitivity_summary.csv"))
saveRDS(sensitivity, here::here("data", "processed", "sensitivity_results.rds"))

message("\n  Sensitivity summary:")
print(sens_summary, n = 20)
}

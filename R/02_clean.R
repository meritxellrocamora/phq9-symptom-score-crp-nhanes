source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in   <- here::here("data", "raw", "download_manifest.csv")
.self <- here::here("R", "02_clean.R")
.out  <- here::here("data", "processed", "analytic_dataset.rds")
if (file.exists(.out) &&
    file.mtime(.out) > file.mtime(.in) &&
    file.mtime(.out) > file.mtime(.self)) {
  message("  Up to date, skipping")
} else {

read_tbl <- function(tbl, id) {
  f1 <- here::here("data", "raw", paste0(tbl, "_", id, ".rds"))
  f2 <- here::here("data", "raw", paste0(id, "_", tbl, ".rds"))
  if (file.exists(f1)) readRDS(f1)
  else if (file.exists(f2)) readRDS(f2)
  else NULL
}

message("  Cleaning and merging across cycles...")

cycle_results <- pmap(cycles, function(id, label, duration, wt_blood, ...) {
  message("Processing ", label)

  demo <- read_tbl("DEMO", id)
  dpq  <- read_tbl("DPQ", id)
  if (is.null(demo) || is.null(dpq)) return(NULL)

  fc <- list(cycle = label, n_demo = nrow(demo))

  # WTPH2YR (phlebotomy weight, 2021-2023) lives in HSCRP_L, NOT in DEMO_L.
  if (wt_blood == "WTPH2YR" && !(wt_blood %in% names(demo))) {
    lab_with_wtph <- read_tbl("HSCRP", id)
    if (!is.null(lab_with_wtph) && "WTPH2YR" %in% names(lab_with_wtph)) {
      demo <- demo |> left_join(
        lab_with_wtph |> select(SEQN, WTPH2YR), by = "SEQN")
      message("  WTPH2YR joined from HSCRP (", id, ")")
    }
  }

  if (!(wt_blood %in% names(demo)))
    stop("Cycle ", label, ": weight variable ", wt_blood, " not found")

  # Marital status column (DMDMARTL in older, DMDMARTZ in newer)
  mart_col <- if ("DMDMARTZ" %in% names(demo)) "DMDMARTZ"
              else if ("DMDMARTL" %in% names(demo)) "DMDMARTL"
              else NULL

  demo_cols <- c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "DMDEDUC2",
                 "INDFMPIR", "SDMVPSU", "SDMVSTRA", "RIDEXPRG", wt_blood)
  if (!is.null(mart_col)) demo_cols <- c(demo_cols, mart_col)

  df <- demo |>
    select(any_of(demo_cols)) |>
    rename(wt_blood_orig = !!sym(wt_blood)) |>
    filter(RIDAGEYR >= 20, wt_blood_orig > 0)

  if (!is.null(mart_col))
    df <- df |>
      rename(marital_raw = !!sym(mart_col)) |>
      mutate(marital_src = mart_col)

  fc$n_under_18 <- sum(demo$RIDAGEYR < 18, na.rm = TRUE)
  fc$n_age_18_19 <- sum(demo$RIDAGEYR %in% 18:19)
  fc$n_no_phleb_wt <- if (wt_blood == "WTPH2YR")
    sum(demo$RIDAGEYR >= 20 & demo$WTMEC2YR > 0 &
          coalesce(demo$WTPH2YR, 0) <= 0, na.rm = TRUE) else 0L
  fc$n_weighted <- nrow(df)
  # Pregnant participants (RIDEXPRG = 1, women aged 20-44) excluded: pregnancy
  # raises CRP and changes sleep, appetite, weight and lipids
  fc$n_pregnant <- sum(df$RIDEXPRG == 1, na.rm = TRUE)
  df <- df |> filter(!(RIDEXPRG %in% 1)) |> select(-RIDEXPRG)
  fc$n_eligible <- nrow(df)

  # --- PHQ-9 + AES ---
  dpq <- dpq |>
    select(SEQN, all_of(phq9_vars)) |>
    mutate(
      across(all_of(phq9_vars), ~ case_when(
        .x %in% 0:3 ~ as.integer(.x),
        TRUE ~ NA_integer_
      )),
      phq9_total = rowSums(pick(all_of(phq9_vars)), na.rm = FALSE),
      aes_score  = rowSums(pick(all_of(aes_items)), na.rm = FALSE),
      depression = as.integer(phq9_total >= depression_cutoff)
    )

  if (mean(!is.na(dpq$phq9_total)) < 0.3)
    warning("Cycle ", label, ": PHQ-9 completion unusually low (",
            round(100 * mean(!is.na(dpq$phq9_total)), 1), "%)")

  df <- left_join(df, dpq, by = "SEQN")
  fc$n_phq9_complete <- sum(!is.na(df$phq9_total))

  # --- CRP with bridging ---
  hscrp <- read_tbl("HSCRP", id)
  if (is.null(hscrp) || !("LBXHSCRP" %in% names(hscrp)) || nrow(hscrp) == 0) {
    hscrp <- read_tbl("CRP", id)
    if (!is.null(hscrp) && "LBXCRP" %in% names(hscrp))
      hscrp <- hscrp |> rename(LBXHSCRP = LBXCRP)
  }
  if (!is.null(hscrp) && "LBXHSCRP" %in% names(hscrp)) {
    hscrp <- hscrp |> select(SEQN, crp_native = LBXHSCRP)
    if (id == "L") {
      hscrp <- hscrp |>
        mutate(crp = crp_xwalk_int + crp_xwalk_slope * crp_native)
    } else if (id == "I") {
      hscrp <- hscrp |>
        mutate(crp = if_else(crp_native <= crp_xwalk_I_max,
                             crp_xwalk_I_int + crp_xwalk_I_slope * crp_native,
                             crp_native))
    } else {
      hscrp <- hscrp |> mutate(crp = crp_native)
    }
    df <- left_join(df, hscrp, by = "SEQN")
    crp_med <- median(df$crp, na.rm = TRUE)
    if (!is.na(crp_med) && crp_med > 10)
      stop("Cycle ", label, ": CRP median ", round(crp_med, 2),
           " mg/L implausibly high")
  }

  if ("crp" %in% names(df)) {
    fc$n_crp_acute <- sum(df$crp > crp_acute_cutoff, na.rm = TRUE)
    df <- df |> mutate(crp_raw = crp,
                       crp = if_else(crp > crp_acute_cutoff, NA_real_, crp))
    fc$n_crp_available <- sum(!is.na(df$crp))
  }

  # --- CBC: derive NLR, PLR, SII ---
  cbc <- read_tbl("CBC", id)
  if (!is.null(cbc) &&
      all(c("LBXWBCSI", "LBXLYPCT", "LBXNEPCT", "LBXPLTSI") %in% names(cbc))) {
    cbc_clean <- cbc |>
      select(SEQN, LBXWBCSI, LBXLYPCT, LBXNEPCT, LBXPLTSI) |>
      mutate(
        lymphocyte_abs = LBXWBCSI * LBXLYPCT / 100,
        neutrophil_abs = LBXWBCSI * LBXNEPCT / 100,
        nlr = if_else(lymphocyte_abs > 0,
                      neutrophil_abs / lymphocyte_abs, NA_real_),
        plr = if_else(lymphocyte_abs > 0,
                      LBXPLTSI / lymphocyte_abs, NA_real_),
        sii = if_else(lymphocyte_abs > 0,
                      LBXPLTSI * neutrophil_abs / lymphocyte_abs, NA_real_),
        # Range checks
        nlr = if_else(nlr > 50, NA_real_, nlr),
        plr = if_else(plr > 1000, NA_real_, plr),
        sii = if_else(sii > 5000, NA_real_, sii)
      ) |>
      select(SEQN, wbc = LBXWBCSI, nlr, plr, sii)
    df <- left_join(df, cbc_clean, by = "SEQN")
    fc$n_cbc_available <- sum(!is.na(df$nlr))
  }

  # --- BMI + waist circumference ---
  bmx <- read_tbl("BMX", id)
  if (!is.null(bmx))
    df <- left_join(df, bmx |>
      select(SEQN, any_of(c("BMXBMI", "BMXWAIST"))) |>
      rename(any_of(c(bmi = "BMXBMI", waist = "BMXWAIST"))), by = "SEQN")

  # --- Metabolic markers ---
  ghb <- read_tbl("GHB", id)
  if (!is.null(ghb) && "LBXGH" %in% names(ghb))
    df <- left_join(df, ghb |> select(SEQN, hba1c = LBXGH), by = "SEQN")

  tchol <- read_tbl("TCHOL", id)
  if (!is.null(tchol) && "LBXTC" %in% names(tchol))
    df <- left_join(df, tchol |> select(SEQN, total_chol = LBXTC), by = "SEQN")

  hdl <- read_tbl("HDL", id)
  if (!is.null(hdl)) {
    hdl_var <- intersect(c("LBDHDD", "LBDHDL"), names(hdl))[1]
    if (!is.na(hdl_var))
      df <- left_join(df, hdl |> select(SEQN, hdl_chol = !!sym(hdl_var)), by = "SEQN")
  }

  trigly <- read_tbl("TRIGLY", id)
  if (!is.null(trigly)) {
    trig_var <- intersect(c("LBXTR", "LBXTLG"), names(trigly))[1]
    if (!is.na(trig_var))
      df <- left_join(df, trigly |>
        select(SEQN, triglycerides = !!sym(trig_var), WTSAF2YR), by = "SEQN")
  }

  # --- Blood pressure harmonization ---
  bpx <- read_tbl("BPX", id)
  if (is.null(bpx)) bpx <- read_tbl("BPXO", id)
  if (!is.null(bpx)) {
    sbp_cols <- intersect(c("BPXSY1","BPXSY2","BPXSY3","BPXSY4",
                             "BPXOSY1","BPXOSY2","BPXOSY3"), names(bpx))
    dbp_cols <- intersect(c("BPXDI1","BPXDI2","BPXDI3","BPXDI4",
                             "BPXODI1","BPXODI2","BPXODI3"), names(bpx))
    if (any(grepl("^BPXSY\\d", sbp_cols)) && any(grepl("^BPXOSY", sbp_cols))) {
      sbp_cols <- grep("^BPXOSY", sbp_cols, value = TRUE)
      dbp_cols <- grep("^BPXODI", dbp_cols, value = TRUE)
    }
    first_sbp <- grep("^BPX(O?SY)1$", sbp_cols, value = TRUE)
    first_dbp <- grep("^BPX(O?DI)1$", dbp_cols, value = TRUE)
    sbp_avg <- setdiff(sbp_cols, first_sbp)
    dbp_avg <- setdiff(dbp_cols, first_dbp)
    if (length(sbp_avg) == 0) sbp_avg <- sbp_cols
    if (length(dbp_avg) == 0) dbp_avg <- dbp_cols
    bp_mat <- as.data.frame(bpx[c("SEQN", sbp_cols, dbp_cols)])
    nan_to_na <- function(x) replace(x, is.nan(x), NA_real_)
    bp_out <- data.frame(SEQN = bp_mat$SEQN)
    if (length(sbp_avg) > 0)
      bp_out$sbp <- nan_to_na(rowMeans(bp_mat[sbp_avg], na.rm = TRUE))
    if (length(dbp_avg) > 0)
      bp_out$dbp <- nan_to_na(rowMeans(bp_mat[dbp_avg], na.rm = TRUE))
    df <- left_join(df, bp_out, by = "SEQN")
  }

  # --- Covariates (all numeric codes, translate = FALSE) ---

  # Smoking
  smq <- read_tbl("SMQ", id)
  if (!is.null(smq) && all(c("SMQ020", "SMQ040") %in% names(smq))) {
    df <- left_join(df, smq |>
      mutate(smoking = factor(case_when(
        SMQ040 %in% c(1, 2) ~ "Current",
        SMQ020 == 1 & SMQ040 == 3 ~ "Former",
        SMQ020 == 2 ~ "Never"
      ), c("Never", "Former", "Current"))) |>
      select(SEQN, smoking), by = "SEQN")
  }

  # Alcohol: any drink in the past 12 months, harmonised across cycles.
  # 2015-2016: ALQ120Q (days per week/month/year; 0 = never) is asked of
  #   ALQ101 = 1 (>= 12 drinks in any one year) or ALQ110 = 1 (>= 12 drinks
  #   in lifetime); ALQ101 = 2 and ALQ110 = 2 (< 12 lifetime drinks) = no
  #   current use. 2017-2018 and 2021-2023: ALQ121 (0 = never in the last
  #   year, 1-10 = drank) is asked of ALQ111 = 1; ALQ111 = 2 = no current use.
  alq <- read_tbl("ALQ", id)
  if (!is.null(alq)) {
    if (all(c("ALQ101", "ALQ110", "ALQ120Q") %in% names(alq))) {
      alq <- alq |> mutate(alcohol = case_when(
        ALQ120Q %in% 1:366 ~ 1L,
        ALQ120Q == 0 ~ 0L,
        ALQ101 == 2 & ALQ110 == 2 ~ 0L))
    } else if (all(c("ALQ111", "ALQ121") %in% names(alq))) {
      alq <- alq |> mutate(alcohol = case_when(
        ALQ121 %in% 1:10 ~ 1L,
        ALQ121 == 0 ~ 0L,
        ALQ111 == 2 ~ 0L))
    } else {
      stop("Cycle ", label, ": alcohol items not recognised")
    }
    df <- left_join(df, alq |> select(SEQN, alcohol), by = "SEQN")
  }

  # Physical activity: leisure-time MVPA (min/week, vigorous counted double),
  # harmonised across cycles. Work and transport domains were not collected
  # in 2021-2023, so only the recreational / leisure-time domain is used.
  # 2015-2018 (GPAQ): PAQ650/PAQ655/PAD660 vigorous recreational and
  #   PAQ665/PAQ670/PAD675 moderate recreational (days/week x min/day).
  # 2021-2023: PAD810Q/PAD810U/PAD820 vigorous LTPA and PAD790Q/PAD790U/
  #   PAD800 moderate LTPA (frequency per day/week/month/year x min/occasion;
  #   codebook: PAD790 = moderate, PAD810 = vigorous).
  paq <- read_tbl("PAQ", id)
  if (!is.null(paq)) {
    na_codes <- \(codes) \(x) if_else(x %in% codes, NA_real_, as.numeric(x))
    gpaq_vars <- c("PAQ650", "PAQ655", "PAD660", "PAQ665", "PAQ670", "PAD675")
    ltpa_vars <- c("PAD790Q", "PAD790U", "PAD800", "PAD810Q", "PAD810U",
                   "PAD820")
    if (all(gpaq_vars %in% names(paq))) {
      pa_df <- paq |>
        select(SEQN, all_of(gpaq_vars)) |>
        mutate(across(c(PAQ650, PAQ665), na_codes(c(7, 9))),
               across(c(PAQ655, PAQ670), na_codes(c(77, 99))),
               across(c(PAD660, PAD675), na_codes(c(7777, 9999))),
               vig_rec = case_when(PAQ650 == 2 ~ 0,
                                   PAQ650 == 1 ~ PAQ655 * PAD660),
               mod_rec = case_when(PAQ665 == 2 ~ 0,
                                   PAQ665 == 1 ~ PAQ670 * PAD675))
    } else if (all(ltpa_vars %in% names(paq))) {
      per_week <- c(D = 7, W = 1, M = 12 / 52, Y = 1 / 52)
      pa_df <- paq |>
        select(SEQN, all_of(ltpa_vars)) |>
        mutate(across(c(PAD790Q, PAD800, PAD810Q, PAD820),
                      na_codes(c(7777, 9999))),
               mod_rec = case_when(
                 PAD790Q == 0 ~ 0,
                 TRUE ~ PAD790Q * unname(per_week[PAD790U]) * PAD800),
               vig_rec = case_when(
                 PAD810Q == 0 ~ 0,
                 TRUE ~ PAD810Q * unname(per_week[PAD810U]) * PAD820))
    } else {
      stop("Cycle ", label, ": physical activity items not recognised")
    }
    pa_df <- pa_df |>
      mutate(pa_leisure_min = 2 * vig_rec + mod_rec,
             pa_meets_guidelines = as.integer(pa_leisure_min >= 150)) |>
      select(SEQN, pa_leisure_min, pa_meets_guidelines)
    df <- left_join(df, pa_df, by = "SEQN")
  }

  # Diabetes
  diq <- read_tbl("DIQ", id)
  if (!is.null(diq) && "DIQ010" %in% names(diq)) {
    df <- left_join(df, diq |>
      mutate(diabetes = case_when(DIQ010 == 1 ~ 1L, DIQ010 %in% c(2, 3) ~ 0L)) |>
      select(SEQN, diabetes), by = "SEQN")
  }

  # Health insurance
  hiq <- read_tbl("HIQ", id)
  if (!is.null(hiq) && "HIQ011" %in% names(hiq)) {
    df <- left_join(df, hiq |>
      mutate(insured = case_when(HIQ011 == 1 ~ 1L, HIQ011 == 2 ~ 0L)) |>
      select(SEQN, insured), by = "SEQN")
  }

  fc$n_final <- nrow(df)
  list(data = df |> mutate(cycle = label, cycle_duration = duration),
       flow = as_tibble(fc))
})

# ---- Stack cycles ----
cycle_results <- compact(cycle_results)
if (length(cycle_results) < 2)
  stop("Only ", length(cycle_results), " cycles loaded; check data/")
message("  ", length(cycle_results), " of ", nrow(cycles), " cycles loaded. Stacking...")

pooled <- bind_rows(map(cycle_results, "data")) |>
  mutate(
    wt_mec = wt_blood_orig / n_cycles,
    wt_fasting = coalesce(WTSAF2YR, 0) / n_cycles,
    sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
    race_ethnicity = factor(case_when(
      RIDRETH3 == 3 ~ "NH White",
      RIDRETH3 == 4 ~ "NH Black",
      RIDRETH3 == 6 ~ "NH Asian",
      RIDRETH3 == 1 ~ "Mexican American",
      RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 7 ~ "Other/Multi-Racial"
    ), c("NH White", "NH Black", "NH Asian",
         "Mexican American", "Other Hispanic", "Other/Multi-Racial")),
    education = factor(case_when(
      DMDEDUC2 %in% c(1, 2) ~ "Less than HS",
      DMDEDUC2 == 3 ~ "HS/GED",
      DMDEDUC2 == 4 ~ "Some college",
      DMDEDUC2 == 5 ~ "College+"
    ), c("Less than HS", "HS/GED", "Some college", "College+")),
    age_group = factor(case_when(
      RIDAGEYR <= 39 ~ "20-39",
      RIDAGEYR <= 59 ~ "40-59",
      TRUE ~ "60+"
    ), c("20-39", "40-59", "60+")),
    log_crp = log(crp + crp_log_offset)
  )

# Marital status, recoded by source column.
# DMDMARTL (2015-2018): 1=Married, 2=Widowed, 3=Divorced, 4=Separated,
#   5=Never married, 6=Living with partner, 77=Refused, 99=DK
# DMDMARTZ (2021-2023): 1=Married/Living with partner,
#   2=Widowed/Divorced/Separated, 3=Never married, 77=Refused, 99=DK
if ("marital_raw" %in% names(pooled)) {
  pooled <- pooled |> mutate(
    marital = factor(case_when(
      marital_src == "DMDMARTL" & marital_raw %in% c(1, 6) ~ "Married/Partner",
      marital_src == "DMDMARTL" & marital_raw %in% c(2, 3, 4) ~
        "Widowed/Divorced/Separated",
      marital_src == "DMDMARTL" & marital_raw == 5 ~ "Never married",
      marital_src == "DMDMARTZ" & marital_raw == 1 ~ "Married/Partner",
      marital_src == "DMDMARTZ" & marital_raw == 2 ~
        "Widowed/Divorced/Separated",
      marital_src == "DMDMARTZ" & marital_raw == 3 ~ "Never married"
    ), c("Married/Partner", "Widowed/Divorced/Separated", "Never married"))
  )
}

# --- Model covariate guards ---
stopifnot(all(model_covars %in% names(pooled)))
covar_cycle_na <- pooled |>
  group_by(cycle) |>
  summarise(across(all_of(model_covars_pa), ~ all(is.na(.x))), .groups = "drop") |>
  pivot_longer(-cycle, names_to = "covariate", values_to = "all_missing") |>
  filter(all_missing)
if (nrow(covar_cycle_na) > 0)
  stop("Model covariate 100% missing within cycle: ",
       paste0(covar_cycle_na$covariate, " (", covar_cycle_na$cycle, ")",
              collapse = "; "), call. = FALSE)

# Per-cycle consistency: weighted prevalence of each categorical covariate
# level should not differ by more than 10 points between cycles
covar_drift <- map_dfr(
  c("sex", "race_ethnicity", "education", "smoking", "alcohol",
    "pa_meets_guidelines"),
  \(v) {
    pooled |>
      filter(!is.na(.data[[v]])) |>
      mutate(level = as.character(.data[[v]])) |>
      group_by(cycle, level) |>
      summarise(w = sum(wt_mec), .groups = "drop_last") |>
      mutate(prop = w / sum(w)) |>
      group_by(level) |>
      summarise(spread = max(prop) - min(prop), .groups = "drop") |>
      mutate(covariate = v)
  }) |>
  filter(spread > 0.10)
drift_known <- covar_drift |> filter(covariate %in% known_cycle_discontinuity)
drift_new   <- covar_drift |> filter(!covariate %in% known_cycle_discontinuity)
if (nrow(drift_known) > 0)
  message("  Known instrument discontinuity between cycles: ",
          paste0(drift_known$covariate, "=", drift_known$level, " (",
                 round(100 * drift_known$spread, 1), " pts)", collapse = "; "))
if (nrow(drift_new) > 0)
  warning("Covariate prevalence differs by > 10 points between cycles: ",
          paste0(drift_new$covariate, "=", drift_new$level, " (",
                 round(100 * drift_new$spread, 1), " pts)",
                 collapse = "; "), call. = FALSE)

# --- IMD classification ---
pooled <- pooled |> mutate(
  dep_subtype = factor(case_when(
    is.na(phq9_total) | phq9_total < depression_cutoff ~ "No depression",
    aes_score >= aes_threshold & crp > crp_cutoff_primary ~ "IMD",
    TRUE ~ "Non-IMD depression"
  ), c("No depression", "Non-IMD depression", "IMD")),
  dep_subtype_strict = factor(case_when(
    is.na(phq9_total) | phq9_total < depression_cutoff ~ "No depression",
    aes_score >= aes_threshold & crp > crp_cutoff_strict ~ "IMD",
    TRUE ~ "Non-IMD depression"
  ), c("No depression", "Non-IMD depression", "IMD")),
  imd_binary = as.integer(dep_subtype == "IMD")
)

# --- Validation guards ---
stopifnot(!any(duplicated(pooled$SEQN)))
stopifnot(all(is.na(pooled$phq9_total) | pooled$phq9_total %in% 0:27))
stopifnot(all(is.na(pooled$aes_score) | pooled$aes_score %in% 0:12))

race_pct <- prop.table(table(pooled$race_ethnicity, useNA = "no"))
if (any(race_pct > 0.80))
  warning("Race distribution skewed: ",
          paste(names(race_pct[race_pct > 0.80]), collapse = ", "), " > 80%")

# --- Exclusion flowchart (unweighted counts) ---
message("\n--- Exclusion Flowchart ---")
flow <- list()
flow$n_adults <- sum(map_dbl(cycle_results,
                             \(x) x$flow$n_demo - x$flow$n_under_18))
flow$n_age_18_19 <- sum(map_dbl(cycle_results, \(x) x$flow$n_age_18_19))
flow$n_no_phleb_wt <- sum(map_dbl(cycle_results,
                                  \(x) x$flow$n_no_phleb_wt))
flow$n_weighted <- sum(map_dbl(cycle_results, \(x) x$flow$n_weighted))
flow$n_pregnant <- sum(map_dbl(cycle_results, \(x) x$flow$n_pregnant))
flow$n_total <- nrow(pooled)
message("Adults 18+ interviewed: ", flow$n_adults)
message("Aged 18-19 excluded (education not asked): ", flow$n_age_18_19)
message("Adults 20+ with weight > 0: ", flow$n_weighted,
        "; pregnant excluded: ", flow$n_pregnant)
message("Eligible adults (20+, non-pregnant, weight > 0): ", flow$n_total)

pooled_analytic <- pooled |> filter(!is.na(phq9_total))
flow$n_phq9 <- nrow(pooled_analytic)
message("After excluding missing PHQ-9: ", flow$n_phq9,
        " (excluded ", flow$n_total - flow$n_phq9, ")")

flow$n_crp_missing <- sum(is.na(pooled_analytic$crp_raw))
flow$n_crp_acute <- sum(pooled_analytic$crp_raw > crp_acute_cutoff,
                        na.rm = TRUE)
pooled_analytic <- pooled_analytic |> filter(!is.na(crp))
flow$n_crp <- nrow(pooled_analytic)
message("After excluding missing/acute CRP: ", flow$n_crp,
        " (excluded ", flow$n_phq9 - flow$n_crp, ")")

key_covars <- c("RIDAGEYR", "sex", "race_ethnicity", "education", "INDFMPIR")
pooled_analytic <- pooled_analytic |>
  filter(if_all(all_of(key_covars), ~ !is.na(.x)))
flow$n_covars <- nrow(pooled_analytic)
message("After excluding missing key covariates: ", flow$n_covars,
        " (excluded ", flow$n_crp - flow$n_covars, ")")

flow_df <- as_tibble(flow)
dir.create(here::here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
write_csv(flow_df, here::here("output", "tables", "exclusion_flowchart.csv"))

# --- Missing data characterization ---
message("\n--- Missing Data ---")
pooled_excluded <- pooled |> filter(!SEQN %in% pooled_analytic$SEQN)

message("Per-variable missingness (full sample):")
miss_vars <- c("phq9_total", "crp", "bmi", "hba1c", "total_chol",
               "hdl_chol", "triglycerides", "nlr", "smoking")
for (v in intersect(miss_vars, names(pooled))) {
  pct <- round(100 * mean(is.na(pooled[[v]])), 1)
  message("  ", v, ": ", pct, "%")
}

if (nrow(pooled_excluded) > 0) {
  message("\nIncluded vs excluded demographics:")
  pooled$included <- pooled$SEQN %in% pooled_analytic$SEQN

  incl_excl <- tibble(
    variable = c("N", "Age (mean)", "PIR (mean)", "Female (%)",
                 names(table(pooled_analytic$race_ethnicity)) |>
                   paste0(" (%)"),
                 names(table(pooled_analytic$education)) |>
                   paste0(" (%)"),
                 "Current smoker (%)", "Depression (%)"),
    included = c(
      sum(pooled$included),
      round(mean(pooled_analytic$RIDAGEYR, na.rm = TRUE), 1),
      round(mean(pooled_analytic$INDFMPIR, na.rm = TRUE), 2),
      round(100 * mean(pooled_analytic$sex == "Female", na.rm = TRUE), 1),
      round(100 * prop.table(table(pooled_analytic$race_ethnicity)), 1),
      round(100 * prop.table(table(pooled_analytic$education)), 1),
      round(100 * mean(pooled_analytic$smoking == "Current", na.rm = TRUE), 1),
      round(100 * mean(pooled_analytic$depression == 1, na.rm = TRUE), 1)
    ),
    excluded = c(
      sum(!pooled$included),
      round(mean(pooled_excluded$RIDAGEYR, na.rm = TRUE), 1),
      round(mean(pooled_excluded$INDFMPIR, na.rm = TRUE), 2),
      round(100 * mean(pooled_excluded$sex == "Female", na.rm = TRUE), 1),
      if (sum(!is.na(pooled_excluded$race_ethnicity)) > 0)
        round(100 * prop.table(table(pooled_excluded$race_ethnicity)), 1)
      else rep(NA_real_, nlevels(pooled_analytic$race_ethnicity)),
      if (sum(!is.na(pooled_excluded$education)) > 0)
        round(100 * prop.table(table(pooled_excluded$education)), 1)
      else rep(NA_real_, nlevels(pooled_analytic$education)),
      round(100 * mean(pooled_excluded$smoking == "Current", na.rm = TRUE), 1),
      round(100 * mean(pooled_excluded$depression == 1, na.rm = TRUE), 1)
    )
  )
  write_csv(incl_excl, here::here("output", "tables", "included_vs_excluded.csv"))
  message("  Included vs excluded comparison saved to output/tables/included_vs_excluded.csv")

  pooled$included <- NULL
}

# --- Save ---
dir.create(here::here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(pooled_analytic, here::here("data", "processed", "analytic_dataset.rds"))
saveRDS(pooled, here::here("data", "processed", "full_dataset.rds"))

flow_by_cycle <- bind_rows(map(cycle_results, "flow"))
saveRDS(flow_by_cycle, here::here("data", "processed", "flow_data.rds"))

# --- Survey design: build once, reuse everywhere ---
design <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~wt_mec, nest = TRUE,
  data = pooled_analytic
)
design_dep <- subset(design, depression == 1)
saveRDS(design, here::here("data", "processed", "survey_design.rds"))
saveRDS(design_dep, here::here("data", "processed", "survey_design_dep.rds"))

message("\n  Analytic N = ", nrow(pooled_analytic))
message("  Depressed (PHQ-9 >= 10): ", sum(pooled_analytic$depression == 1, na.rm = TRUE))
message("  IMD (primary): ", sum(pooled_analytic$dep_subtype == "IMD", na.rm = TRUE))
message("  IMD (strict): ", sum(pooled_analytic$dep_subtype_strict == "IMD", na.rm = TRUE))
}

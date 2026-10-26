source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in   <- here::here("data", "processed",
                    c("full_dataset.rds", "analytic_dataset.rds",
                      "analysis_results.rds"))
.self <- here::here("R", "06_sensitivity_ipcw_mi.R")
.out  <- here::here("data", "processed", "ipcw_mi_results.rds")
if (file.exists(.out) &&
    file.mtime(.out) > max(file.mtime(.in)) &&
    file.mtime(.out) > file.mtime(.self)) {
  message("  Up to date, skipping")
} else {

set.seed(12345)
full_df <- readRDS(here::here("data", "processed", "full_dataset.rds"))
analytic <- readRDS(here::here("data", "processed", "analytic_dataset.rds"))

message("  Full dataset N (eligible adults 20+): ", nrow(full_df))
message("  Analytic N (complete-case): ", nrow(analytic))
message("  Excluded N: ", nrow(full_df) - nrow(analytic))

# ============================================================
# 1. IPCW: predict probability of inclusion
# ============================================================
message("\n  [1/4] IPCW: modelling P(included | pre-exclusion covariates)...")

# CRP above the acute-inflammation cut-off is a design exclusion, not
# nonresponse, so those adults are outside the inclusion model.
ipcw_df <- full_df |>
  filter(is.na(crp_raw) | crp_raw <= crp_acute_cutoff) |>
  mutate(included = as.integer(SEQN %in% analytic$SEQN))
message("  Inclusion-model denominator (CRP <= ", crp_acute_cutoff,
        " mg/L or missing): ", nrow(ipcw_df))

# Predictors: fully observed pre-exclusion variables only (age, sex,
# race/ethnicity, cycle). Education and PIR are themselves inclusion
# criteria, so their missingness indicators would perfectly predict
# exclusion. Assumed mechanism: inclusion is MAR given these four
# variables; selection related to SES or to depression itself is not
# corrected by this model (stated in Limitations).
ipcw_df <- ipcw_df |>
  mutate(cycle_f = factor(cycle))

ipcw_predictors <- c("RIDAGEYR", "sex", "race_ethnicity", "cycle_f")

ipcw_design_full <- svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~wt_mec, nest = TRUE, data = ipcw_df
)

ipcw_formula <- reformulate(ipcw_predictors, "included")
ipcw_mod <- svyglm(ipcw_formula, design = ipcw_design_full,
                   family = quasibinomial())
ipcw_df$p_included <- predict(ipcw_mod, type = "response",
                              newdata = ipcw_df) |> as.numeric()

# Stabilized IPCW: numerator = survey-weighted marginal P(included),
# denominator = P(included | X)
p_marg <- as.numeric(coef(svymean(~included, ipcw_design_full)))
ipcw_df$ipcw <- with(ipcw_df, p_marg / p_included)

# Diagnostic: truncate at 1st/99th percentile to limit influence of extreme weights
q01 <- quantile(ipcw_df$ipcw, 0.01, na.rm = TRUE)
q99 <- quantile(ipcw_df$ipcw, 0.99, na.rm = TRUE)
ipcw_df$ipcw_trunc <- pmin(pmax(ipcw_df$ipcw, q01), q99)

message(sprintf("  Marginal P(included) = %.3f", p_marg))
message(sprintf("  IPCW range: [%.3f, %.3f], median %.3f",
                min(ipcw_df$ipcw, na.rm = TRUE),
                max(ipcw_df$ipcw, na.rm = TRUE),
                median(ipcw_df$ipcw, na.rm = TRUE)))
message(sprintf("  After 1%%/99%% truncation: [%.3f, %.3f]", q01, q99))

# Diagnostic: c-statistic of inclusion model
roc_auc <- tryCatch({
  pr <- ipcw_df$p_included
  y  <- ipcw_df$included
  rk <- rank(pr)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(rk[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}, error = function(e) NA_real_)
message(sprintf("  IPCW model c-statistic: %.3f", roc_auc))

# Merge IPCW back onto analytic data
analytic_ipcw <- analytic |>
  left_join(ipcw_df |> select(SEQN, ipcw, ipcw_trunc), by = "SEQN")

# ============================================================
# 2. Multiple imputation (mice) on analytic sample
# ============================================================
message("\n  [2/4] Multiple imputation with mice (20 imputations)...")

# Variables to impute: any covariate with missingness in analytic set.
# Outcomes (imd_binary, depression, dep_subtype) are deterministic functions
# of phq9_total/aes_score/crp and are NOT imputed; pass them as "auxiliary"
# fully-observed predictors (mice will treat as observed).

mice_input <- analytic_ipcw |>
  select(SEQN, RIDAGEYR, sex, race_ethnicity, education, INDFMPIR,
         smoking, alcohol, pa_meets_guidelines,
         phq9_total, aes_score, crp, log_crp, bmi, hba1c,
         depression, imd_binary, dep_subtype,
         age_group, SDMVPSU, SDMVSTRA, wt_mec, ipcw_trunc, cycle)

# Missingness in mice input (within analytic set)
miss_pct <- map_dbl(mice_input, ~ 100 * mean(is.na(.x)))
miss_pct <- miss_pct[miss_pct > 0]
message("  Covariates requiring imputation (analytic set):")
for (v in names(miss_pct)) {
  message(sprintf("    %-25s %5.1f%%", v, miss_pct[v]))
}

# Predictor matrix: survey cycle and the analysis weight enter every
# imputation model (design-informed imputation); SEQN, strata, PSU and the
# IPCW weight do not. Variables without missing values are not imputed.
pred_init <- mice::quickpred(mice_input, mincor = 0.1)
no_pred <- c("SEQN", "SDMVPSU", "SDMVSTRA", "ipcw_trunc")
pred_init[no_pred, ] <- 0
pred_init[, no_pred] <- 0
pred_init[, c("cycle", "wt_mec")] <- 1
# depression is a deterministic function of dep_subtype / imd_binary
pred_init[, "depression"] <- 0
diag(pred_init) <- 0
mice_input <- mice_input |> mutate(cycle = factor(cycle))

# Method assignment: pmm for continuous, logreg/polyreg for categorical;
# variables without missing values are not imputed
meth_init <- mice::make.method(mice_input)
for (v in names(mice_input)) {
  if (!anyNA(mice_input[[v]])) meth_init[v] <- ""
}

message(sprintf("  Variables actively imputed: %d",
                sum(meth_init != "")))

# Run mice
m_imp <- 20L
mice_fit <- mice::mice(mice_input,
                      m = m_imp,
                      method = meth_init,
                      predictorMatrix = pred_init,
                      maxit = 10L,
                      printFlag = FALSE,
                      seed = 12345)

# Convergence checks
mice_loggedEvents <- mice_fit$loggedEvents
if (!is.null(mice_loggedEvents) && nrow(mice_loggedEvents) > 0) {
  message("  mice logged events (last 10):")
  print(utils::tail(mice_loggedEvents, 10))
} else {
  message("  mice logged events: none (clean convergence)")
}

# ============================================================
# 3. Parallel survey-weighted analyses across 20 imputations
# ============================================================
message("\n  [3/4] Parallel analysis across 20 imputed datasets...")

plan(multisession, workers = 4L)
on.exit(plan(sequential), add = TRUE)

stopifnot(all(model_covars_pa %in% names(mice_input)))

with_p <- progressr::progressor(steps = m_imp)

fit_one_imputation <- function(i) {
  d <- mice::complete(mice_fit, action = i)
  # IPCW × MEC weight: combined weight for IPCW-MI sensitivity
  d$w_combined <- d$wt_mec * d$ipcw_trunc

  des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                   weights = ~w_combined, nest = TRUE, data = d)
  des_dep <- subset(des, depression == 1)

  # IMD prevalence among depressed
  prev_overall <- svymean(~imd_binary, des_dep, na.rm = TRUE)

  # Predictors model (fully adjusted, same formula as primary)
  pred_formula <- reformulate(model_covars, "imd_binary")
  pred_mod <- svyglm(pred_formula, design = des_dep,
                     family = quasibinomial())

  # Also overall depression prevalence (for completeness)
  dep_overall <- svymean(~depression, des, na.rm = TRUE)

  list(
    prev_imd  = prev_overall,
    prev_dep  = dep_overall,
    df_dep    = survey::degf(des),
    pred_mod  = pred_mod,
    n_model   = nobs(pred_mod),
    df_model  = pred_mod$df.residual,
    df_design = survey::degf(des_dep),
    n_dep     = sum(d$depression == 1, na.rm = TRUE),
    n_imd     = sum(d$imd_binary == 1, na.rm = TRUE)
  )
}

mi_fits <- furrr::future_map(
  seq_len(m_imp),
  function(i) {
    out <- fit_one_imputation(i)
    out
  },
  .options = furrr::furrr_options(seed = FALSE)
)

plan(sequential)

# ============================================================
# 4. Pool via MIcombine
# ============================================================
message("\n  [4/4] Pooling via MIcombine...")

# Small-sample (Barnard-Rubin) adjusted df; dfcom = complete-data model df
pool_rubin <- function(est, v, dfcom) {
  mean_est <- mean(est)
  within <- mean(v)
  between <- var(est)
  total <- within + (1 + 1 / m_imp) * between
  lambda <- max((1 + 1 / m_imp) * between / total, 1e-04)
  df_old <- (m_imp - 1) / lambda^2
  df_obs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
  df <- df_old * df_obs / (df_old + df_obs)
  crit <- qt(0.975, df)
  tibble(estimate = mean_est, se = sqrt(total), df = df,
         CI_lo = max(0, mean_est - crit * sqrt(total)),
         CI_hi = min(1, mean_est + crit * sqrt(total)),
         fmi = (1 + 1 / m_imp) * between / total)
}

# Pool predictor coefficients
beta_list <- map(mi_fits, ~ coef(.x$pred_mod))
vcov_list <- map(mi_fits, ~ vcov(.x$pred_mod))
dfcom_pred <- mi_fits[[1]]$df_model

pooled <- mitools::MIcombine(results = beta_list, variances = vcov_list,
                             df.complete = dfcom_pred)
pooled_se <- sqrt(diag(pooled$variance))
crit <- qt(0.975, df = pooled$df)
# Fraction of missing information as (1 + 1/m) B / T (mitools' missinfo
# adds the small-sample term 2 / (dfcom + 3), which dominates at dfcom = 31)
between_var <- apply(do.call(rbind, beta_list), 2, var)
fmi_lambda <- (1 + 1 / m_imp) * between_var / diag(pooled$variance)
pooled_summary <- tibble(
  term     = names(pooled$coefficients),
  estimate = as.numeric(pooled$coefficients),
  se       = pooled_se,
  OR       = exp(as.numeric(pooled$coefficients)),
  CI_lo    = exp(as.numeric(pooled$coefficients) - crit * pooled_se),
  CI_hi    = exp(as.numeric(pooled$coefficients) + crit * pooled_se),
  fmi      = as.numeric(fmi_lambda),
  df       = pooled$df,
  n_model  = mi_fits[[1]]$n_model,
  df_complete = dfcom_pred
) |>
  mutate(p = 2 * pt(-abs(estimate / se), df = df)) |>
  filter(term != "(Intercept)")

# Pool IMD prevalence among depressed
prev_pooled_imd <- pool_rubin(
  map_dbl(mi_fits, ~ as.numeric(coef(.x$prev_imd))),
  map_dbl(mi_fits, ~ as.numeric(vcov(.x$prev_imd))),
  mi_fits[[1]]$df_design)

message(sprintf("  Pooled IMD prevalence (IPCW + MI): %.1f%% (%.1f-%.1f)",
                100 * prev_pooled_imd$estimate,
                100 * prev_pooled_imd$CI_lo,
                100 * prev_pooled_imd$CI_hi))

# Pool overall depression prevalence
prev_pooled_dep <- pool_rubin(
  map_dbl(mi_fits, ~ as.numeric(coef(.x$prev_dep))),
  map_dbl(mi_fits, ~ as.numeric(vcov(.x$prev_dep))),
  mi_fits[[1]]$df_dep)

# ============================================================
# 5. Comparison to complete-case primary
# ============================================================
results_primary <- readRDS(here::here("data", "processed", "analysis_results.rds"))
cc_prev <- coef(results_primary$prevalence$imd$overall)
cc_prev_ci <- confint(results_primary$prevalence$imd$overall,
                      df = results_primary$prevalence$design_df)
cc_ors <- results_primary$predictors$ors

comparison_prev <- tibble(
  analysis = c("Complete-case (primary)", "IPCW + MI"),
  imd_prev = c(100 * cc_prev, 100 * prev_pooled_imd$estimate),
  ci_lo    = c(100 * cc_prev_ci[1], 100 * prev_pooled_imd$CI_lo),
  ci_hi    = c(100 * cc_prev_ci[2], 100 * prev_pooled_imd$CI_hi),
  se       = c(sqrt(attr(results_primary$prevalence$imd$overall, "var"))[1],
               prev_pooled_imd$se)
)

# Side-by-side predictor OR comparison
comparison_ors <- cc_ors |>
  select(term, OR_cc = OR, CI_lo_cc = CI_lo, CI_hi_cc = CI_hi,
         se_cc = se, p_cc = p, n_cc = n, df_cc = df) |>
  full_join(
    pooled_summary |>
      select(term, OR_mi = OR, CI_lo_mi = CI_lo, CI_hi_mi = CI_hi,
             se_mi = se, p_mi = p, fmi_coef = fmi, n_mi = n_model, df_mi = df),
    by = "term"
  ) |>
  mutate(
    se_ratio = se_mi / se_cc,
    or_pct_change = 100 * (OR_mi - OR_cc) / OR_cc
  )

# ============================================================
# 6. Save outputs
# ============================================================
ipcw_mi_results <- list(
  ipcw = list(
    model = ipcw_mod,
    weights_summary = tibble(
      min = min(ipcw_df$ipcw, na.rm = TRUE),
      q01 = q01,
      median = median(ipcw_df$ipcw, na.rm = TRUE),
      mean = mean(ipcw_df$ipcw, na.rm = TRUE),
      q99 = q99,
      max = max(ipcw_df$ipcw, na.rm = TRUE),
      n_truncated = sum(ipcw_df$ipcw != ipcw_df$ipcw_trunc, na.rm = TRUE)
    ),
    c_statistic = roc_auc,
    p_marginal = p_marg
  ),
  mice = list(
    m = m_imp,
    method = meth_init[meth_init != ""],
    loggedEvents = mice_loggedEvents,
    chain_mean = mice_fit$chainMean,
    chain_var = mice_fit$chainVar
  ),
  pooled = list(
    prev_imd = prev_pooled_imd,
    prev_dep = prev_pooled_dep,
    predictors = pooled_summary
  ),
  comparison = list(
    prevalence = comparison_prev,
    predictors = comparison_ors
  ),
  n = list(
    full = nrow(full_df),
    ipcw_denominator = nrow(ipcw_df),
    analytic = nrow(analytic),
    excluded = nrow(full_df) - nrow(analytic),
    m_imp = m_imp,
    workers = 4L
  )
)

saveRDS(ipcw_mi_results, here::here("data", "processed", "ipcw_mi_results.rds"))

# Write summary CSV for table integration
sens_ipcw_mi <- tibble(
  analysis = "IPCW + MI (m = 20)",
  n_total = nrow(analytic),
  n_depressed = round(mean(map_dbl(mi_fits, "n_dep"))),
  n_imd = round(mean(map_dbl(mi_fits, "n_imd"))),
  imd_prev = round(100 * prev_pooled_imd$estimate, 1),
  imd_ci = sprintf("(%.1f-%.1f)",
                   100 * prev_pooled_imd$CI_lo,
                   100 * prev_pooled_imd$CI_hi),
  df = round(prev_pooled_imd$df, 1)
)
write_csv(sens_ipcw_mi,
          here::here("output", "tables", "sensitivity_ipcw_mi.csv"))

write_csv(comparison_ors,
          here::here("output", "tables", "ipcw_mi_vs_cc_predictors.csv"))

write_csv(tibble(n_denominator = ipcw_mi_results$n$ipcw_denominator,
                 c_statistic = ipcw_mi_results$ipcw$c_statistic,
                 trunc_lo = q01, trunc_hi = q99),
          here::here("output", "tables", "ipcw_inclusion_model.csv"))

message("\n  IPCW + MI results saved to data/processed/ipcw_mi_results.rds")
message("  Summary table: output/tables/sensitivity_ipcw_mi.csv")
message("  Predictor comparison: output/tables/ipcw_mi_vs_cc_predictors.csv")
}

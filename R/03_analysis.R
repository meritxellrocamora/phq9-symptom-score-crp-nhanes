source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in   <- here::here("data", "processed", "analytic_dataset.rds")
.self <- here::here("R", "03_analysis.R")
.out  <- here::here("data", "processed", "analysis_results.rds")
if (file.exists(.out) &&
    file.mtime(.out) > file.mtime(.in) &&
    file.mtime(.out) > file.mtime(.self)) {
  message("  Up to date, skipping")
} else {

df <- readRDS(here::here("data", "processed", "analytic_dataset.rds"))
stopifnot(nrow(df) > 0,
          all(c("depression", "imd_binary", "dep_subtype", "wt_mec",
                "SDMVPSU", "SDMVSTRA") %in% names(df)))
message("  Analytic N = ", nrow(df))

# --- Survey design (built once in 02_clean.R) ---
design     <- readRDS(here::here("data", "processed", "survey_design.rds"))
design_dep <- readRDS(here::here("data", "processed", "survey_design_dep.rds"))

stopifnot(all(c(model_covars_pa, "wt_fasting") %in% names(df)))

# Fasting-subsample design (WTSAF2YR / 3) for triglycerides
design_fast <- fasting_design(df)
design_for  <- \(v) if (v == "triglycerides") design_fast else design

fit <- function(formula_str, des, fam = quasibinomial()) {
  svyglm(as.formula(formula_str), design = des, family = fam)
}

extract_or <- function(mod) {
  s <- summary(mod)$coefficients
  crit <- qt(0.975, mod$df.residual)
  vars <- setdiff(all.vars(formula(mod)), "imd_binary")
  cyc <- paste(sort(unique(df$cycle[df$depression == 1 &
                                     complete.cases(df[vars])])),
               collapse = "; ")
  tibble(
    term     = rownames(s),
    estimate = s[, "Estimate"],
    se       = s[, "Std. Error"],
    OR       = exp(s[, "Estimate"]),
    CI_lo    = exp(s[, "Estimate"] - crit * s[, "Std. Error"]),
    CI_hi    = exp(s[, "Estimate"] + crit * s[, "Std. Error"]),
    p        = s[, "Pr(>|t|)"],
    n        = nobs(mod),
    df       = mod$df.residual,
    cycles   = cyc
  ) |> filter(term != "(Intercept)")
}

# ============================================================
# 1. Depression prevalence
# ============================================================
message("  [1/7] Depression prevalence...")

prev_overall <- svymean(~depression, design, na.rm = TRUE)
prev_by_race <- svyby(~depression, ~race_ethnicity, design, svymean, na.rm = TRUE)
prev_by_sex  <- svyby(~depression, ~sex, design, svymean, na.rm = TRUE)
prev_by_age  <- svyby(~depression, ~age_group, design, svymean, na.rm = TRUE)

message("  Overall depression prevalence: ",
        round(100 * coef(prev_overall), 1), "%")

# ============================================================
# 2. IMD prevalence among depressed
# ============================================================
message("  [2/7] IMD prevalence among depressed...")

imd_overall  <- svymean(~imd_binary, design_dep, na.rm = TRUE)
imd_by_race  <- svyby(~imd_binary, ~race_ethnicity, design_dep, svymean, na.rm = TRUE)
imd_by_sex   <- svyby(~imd_binary, ~sex, design_dep, svymean, na.rm = TRUE)
imd_by_age   <- svyby(~imd_binary, ~age_group, design_dep, svymean, na.rm = TRUE)

# Cell size check
imd_n <- df |> filter(depression == 1) |> count(race_ethnicity, imd_binary)
small_cells <- imd_n |> filter(n < 20)
if (nrow(small_cells) > 0) {
  warning("Small cell sizes (n < 20) in IMD x race:\n",
          paste(capture.output(print(small_cells)), collapse = "\n"))
}

message("  Overall IMD prevalence among depressed: ",
        round(100 * coef(imd_overall), 1), "%")

# Criterion-level proportions among depressed adults: AES >= 6, CRP > 1,
# both, CRP > 1 conditional on AES status, and a Rao-Scott test of
# independence of the two criteria
dep_df <- df |> filter(depression == 1)
design_crit <- update(design_dep,
  aes_met  = as.integer(aes_score >= aes_threshold),
  crp_met  = as.integer(crp > crp_cutoff_primary),
  both_met = as.integer(aes_score >= aes_threshold & crp > crp_cutoff_primary))
crit_means <- svymean(~aes_met + crp_met + both_met, design_crit)
crit_ci <- confint(crit_means, df = survey::degf(design_crit))
crp_by_aes <- svyby(~crp_met, ~aes_met, design_crit, svymean)
crp_by_aes_ci <- confint(crp_by_aes, df = survey::degf(design_crit))
crit_chisq <- svychisq(~aes_met + crp_met, design_crit)
crit_n <- dep_df |>
  summarise(aes_met = sum(aes_score >= aes_threshold),
            crp_met = sum(crp > crp_cutoff_primary),
            both_met = sum(aes_score >= aes_threshold & crp > crp_cutoff_primary),
            crp_met_aes_below = sum(crp > crp_cutoff_primary &
                                      aes_score < aes_threshold),
            n_aes_below = sum(aes_score < aes_threshold))
criterion_overlap <- tibble(
  quantity = c("Symptom score >= 6", "CRP > 1 mg/L", "Both (IMD approximation)",
               "CRP > 1 mg/L among AES >= 6", "CRP > 1 mg/L among AES < 6",
               "Expected both if independent (product of marginals)"),
  n_unweighted = c(crit_n$aes_met, crit_n$crp_met, crit_n$both_met,
                   crit_n$both_met, crit_n$crp_met_aes_below, NA_integer_),
  denominator = c(nrow(dep_df), nrow(dep_df), nrow(dep_df),
                  crit_n$aes_met, crit_n$n_aes_below, nrow(dep_df)),
  pct_weighted = 100 * c(coef(crit_means), coef(crp_by_aes)[c("1", "0")],
                         coef(crit_means)["aes_met"] * coef(crit_means)["crp_met"]),
  ci_lo = 100 * c(crit_ci[, 1], crp_by_aes_ci[c("1", "0"), 1], NA_real_),
  ci_hi = 100 * c(crit_ci[, 2], crp_by_aes_ci[c("1", "0"), 2], NA_real_)
) |>
  mutate(rao_scott_p = c(rep(NA_real_, 5), crit_chisq$p.value),
         across(where(is.double), ~ round(.x, 3)))
# Design-based difference in the CRP proportion (high minus low symptom score)
m_diff <- svyglm(crp_met ~ aes_met, design_crit)
diff_ci <- confint(m_diff, df = m_diff$df.residual)["aes_met", ]
criterion_overlap <- criterion_overlap |>
  add_row(quantity = "Difference, high minus low symptom score (points)",
          n_unweighted = NA_integer_, denominator = NA_integer_,
          pct_weighted = round(100 * coef(m_diff)[["aes_met"]], 1),
          ci_lo = round(100 * diff_ci[[1]], 1), ci_hi = round(100 * diff_ci[[2]], 1),
          rao_scott_p = round(crit_chisq$p.value, 3))
write_csv(criterion_overlap,
          here::here("output", "tables", "table2_criterion_overlap.csv"))

# The same association test for every alternative symptom definition
symptom_defs <- list(
  "Four items, score \u2265 6 (primary)" = \(d) d$aes_score >= 6,
  "Four items, score \u2265 7" = \(d) d$aes_score >= 7,
  "Four items, score \u2265 8" = \(d) d$aes_score >= 8,
  "Fatigue and psychomotor change, score \u2265 3" = \(d) d$DPQ040 + d$DPQ080 >= 3,
  "Sleep, fatigue and appetite (no psychomotor item), score \u2265 5" =
    \(d) d$DPQ030 + d$DPQ040 + d$DPQ050 >= 5,
  "Fatigue item \u2265 2" = \(d) d$DPQ040 >= 2,
  "Fatigue item = 3" = \(d) d$DPQ040 == 3
)
symptom_assoc <- imap_dfr(symptom_defs, \(f, lab) {
  des <- update(design_dep, sym_met = as.integer(f(design_dep$variables)),
                crp_met = as.integer(crp > crp_cutoff_primary))
  by <- svyby(~crp_met, ~sym_met, des, svymean)
  m  <- svyglm(crp_met ~ sym_met, des)
  ci <- confint(m, df = m$df.residual)["sym_met", ]
  chi <- svychisq(~sym_met + crp_met, des)
  dd <- des$variables
  tibble(symptom_definition = lab,
         n_met = sum(dd$sym_met == 1), n_not_met = sum(dd$sym_met == 0),
         pct_met_weighted = round(100 * as.numeric(svymean(~sym_met, des)), 1),
         crp_above1_pct_if_met = round(100 * coef(by)[["1"]], 1),
         crp_above1_pct_if_not_met = round(100 * coef(by)[["0"]], 1),
         difference_points = round(100 * coef(m)[["sym_met"]], 1),
         diff_ci_lo = round(100 * ci[[1]], 1), diff_ci_hi = round(100 * ci[[2]], 1),
         rao_scott_p = round(chi$p.value, 3))
})
write_csv(symptom_assoc,
          here::here("output", "tables", "table2b_symptom_definitions.csv"))

# Four-group metabolic comparison among depressed adults (CRP criterion x
# symptom criterion), with the symptom criterion tested within CRP strata
four_grp <- \(des) update(des,
  crp_met = factor(as.integer(crp > crp_cutoff_primary), c(0, 1),
                   c("CRP \u2264 1 mg/L", "CRP > 1 mg/L")),
  sym_met = factor(as.integer(aes_score >= aes_threshold), c(0, 1),
                   c("score < 6", "score \u2265 6")))
des4 <- four_grp(design_dep)
des4_fast <- four_grp(subset(design_fast, depression == 1))
four_vars <- c("bmi", "waist", "hba1c", "hdl_chol", "triglycerides")
four_group <- map_dfr(four_vars, \(v) {
  d <- if (v == "triglycerides") des4_fast else des4
  means <- svyby(reformulate(v), ~crp_met + sym_met, d, svymean, na.rm = TRUE)
  ns <- d$variables |> filter(!is.na(.data[[v]])) |> count(crp_met, sym_met)
  m_add <- svyglm(reformulate(c("crp_met", "sym_met"), response = v), d)
  m_int <- svyglm(reformulate(c("crp_met", "sym_met", "crp_met:sym_met"),
                              response = v), d)
  tibble(marker = v,
         group = paste(means$crp_met, means$sym_met, sep = ", "),
         n = ns$n[match(paste(means$crp_met, means$sym_met),
                        paste(ns$crp_met, ns$sym_met))],
         mean = as.numeric(coef(means)),
         se = as.numeric(SE(means)),
         p_symptom_adj_crp = regTermTest(m_add, "sym_met")$p[1],
         p_crp_adj_symptom = regTermTest(m_add, "crp_met")$p[1],
         p_interaction = regTermTest(m_int, "crp_met:sym_met")$p[1])
})
write_csv(four_group,
          here::here("output", "tables", "tableS5_four_group_metabolic.csv"))
message("  CRP > 1 among AES >= 6: ",
        round(100 * coef(crp_by_aes)["1"], 1), "% vs AES < 6: ",
        round(100 * coef(crp_by_aes)["0"], 1), "% (Rao-Scott P = ",
        signif(crit_chisq$p.value, 2), ")")

# Subgroup precision (NCHS Data Presentation Standards for Proportions,
# Vital Health Stat 2(175): Korn-Graubard CI, effective n = n / DEff)
dep_df <- df |> filter(depression == 1)
precision_row <- function(des_sub, n_sub, n_imd_sub, stratifier, group) {
  m    <- svymean(~imd_binary, des_sub, deff = TRUE)
  p_l  <- svyciprop(~imd_binary, des_sub, method = "logit")
  p_b  <- svyciprop(~imd_binary, des_sub, method = "beta")
  ci_l <- confint(p_l)
  ci_b <- confint(p_b)
  deff_v <- as.numeric(deff(m))
  eff_n  <- n_sub / deff_v
  width  <- ci_b[2] - ci_b[1]
  rel_width <- width / as.numeric(coef(p_b))
  d_f <- survey::degf(des_sub)
  status <- case_when(
    n_sub < 30 | eff_n < 30 ~ "Does not meet standard (n or effective n below 30)",
    width >= 0.30 ~ "Does not meet standard (interval width 30 points or more)",
    width > 0.05 & rel_width > 1.30 ~ "Does not meet standard (relative width above 130%)",
    d_f < 8 ~ "Flag (fewer than 8 degrees of freedom)",
    TRUE ~ "Meets standard"
  )
  tibble(stratifier = stratifier, group = group,
         n_depressed = n_sub, n_imd = n_imd_sub,
         pct_weighted = 100 * as.numeric(coef(p_l)),
         logit_ci_lo = 100 * ci_l[1], logit_ci_hi = 100 * ci_l[2],
         kg_ci_lo = 100 * ci_b[1], kg_ci_hi = 100 * ci_b[2],
         deff = deff_v, effective_n = eff_n, design_df = d_f,
         ci_width_kg = 100 * width, rel_width_kg = 100 * rel_width,
         nchs_status = status)
}
strat_vars <- c("Race/Ethnicity" = "race_ethnicity", "Sex" = "sex",
                "Age group" = "age_group")
subgroup_precision <- bind_rows(
  precision_row(design_dep, nrow(dep_df), sum(dep_df$imd_binary == 1),
                "Overall", "All depressed adults"),
  map2_dfr(strat_vars, names(strat_vars), \(v, lab) {
    map_dfr(levels(df[[v]]), \(g) {
      des_g <- eval(bquote(subset(design_dep, .(as.name(v)) == .(g))))
      d_g <- dep_df[dep_df[[v]] == g, ]
      precision_row(des_g, nrow(d_g), sum(d_g$imd_binary == 1), lab, g)
    })
  })
)
write_csv(subgroup_precision |> mutate(across(where(is.double), ~ round(.x, 4))),
          here::here("output", "tables", "tableS3_subgroup_precision.csv"))

# ============================================================
# 3. Metabolic profile comparison
# ============================================================
message("  [3/7] Metabolic profile comparison...")

# Triglycerides use the fasting-subsample design (WTSAF2YR / 3, WTSAF2YR > 0).
# Note: CRP comparison is partly tautological (CRP > 1 defines IMD).
# Reported for completeness but should not be interpreted as independent validation.
metabolic_vars <- c("bmi", "waist", "hba1c", "total_chol", "hdl_chol",
                    "triglycerides", "crp", "nlr", "plr", "sii")
metabolic_vars <- intersect(metabolic_vars, names(df))

metabolic_means <- map(metabolic_vars, \(v) {
  tryCatch(
    svyby(reformulate(v), ~dep_subtype, design_for(v),
          svymean, na.rm = TRUE),
    error = function(e) { message("  ", v, ": ", e$message); NULL }
  )
}) |> set_names(metabolic_vars) |> compact()

metabolic_tests <- map(metabolic_vars, \(v) {
  tryCatch(
    svyglm(reformulate("dep_subtype", response = v), design = design_for(v)),
    error = function(e) NULL
  )
}) |> set_names(metabolic_vars) |> compact()

metabolic_info <- tibble(
  marker      = names(metabolic_tests),
  n           = map_int(metabolic_tests, nobs),
  df_residual = map_int(metabolic_tests, \(m) as.integer(m$df.residual))
)

all_p <- map_dbl(metabolic_tests, \(m) {
  tryCatch(regTermTest(m, "dep_subtype")$p[1], error = function(e) NA_real_)
})
# BH family = non-definitional markers only (CRP > 1 mg/L defines IMD)
bh_family <- setdiff(names(all_p), "crp")
bh_adjust <- function(p) {
  out <- p; out[] <- NA_real_
  out[bh_family] <- p.adjust(p[bh_family], method = "BH"); out
}
all_p_adj <- bh_adjust(all_p)

# Pairwise: IMD vs Non-IMD depression (relevel via update, not rebuild)
relevel_imd <- \(d) update(d,
  dep_subtype = relevel(factor(dep_subtype), ref = "Non-IMD depression"))
design_relevel <- relevel_imd(design)
design_fast_relevel <- relevel_imd(design_fast)

pairwise_p <- map_dbl(names(metabolic_tests), \(v) {
  tryCatch({
    m <- svyglm(reformulate("dep_subtype", response = v),
                design = if (v == "triglycerides") design_fast_relevel
                         else design_relevel)
    s <- summary(m)$coefficients
    imd_row <- grep("IMD$", rownames(s))
    if (length(imd_row) == 1) s[imd_row, "Pr(>|t|)"] else NA_real_
  }, error = function(e) NA_real_)
}) |> set_names(names(metabolic_tests))
pairwise_p_adj <- bh_adjust(pairwise_p)

# Prior triglyceride result with examination weights (kept for comparison)
trig_exam <- list(
  means = svyby(~triglycerides, ~dep_subtype, design, svymean, na.rm = TRUE),
  p_omnibus = regTermTest(
    svyglm(triglycerides ~ dep_subtype, design = design),
    "dep_subtype")$p[1],
  p_pairwise = summary(svyglm(triglycerides ~ dep_subtype,
                              design = design_relevel))$coefficients[
    "dep_subtypeIMD", "Pr(>|t|)"]
)

message("  Metabolic markers tested: ", length(metabolic_tests))
message("  Non-definitional markers in BH family: ", length(bh_family))
message("  Significant after BH: ", sum(all_p_adj < 0.05, na.rm = TRUE))
message("  Pairwise IMD vs Non-IMD significant after BH: ",
        sum(pairwise_p_adj < 0.05, na.rm = TRUE), " of ", length(bh_family))

# ============================================================
# 4. Predictors of IMD among depressed
# ============================================================
message("  [4/7] Predictors of IMD...")

demo_covars <- c("RIDAGEYR", "sex", "race_ethnicity")

# Model 1: Crude (each predictor alone)
pred_crude <- map_dfr(model_covars, \(v) {
  m <- svyglm(reformulate(v, "imd_binary"),
              design = design_dep, family = quasibinomial())
  extract_or(m)
}) |>
  mutate(p_adj = p.adjust(p, method = "BH"))

# Model 2: Demographics-adjusted (age, sex, race/ethnicity)
pred_demo_formula <- paste("imd_binary ~", paste(demo_covars, collapse = " + "))
pred_demo <- svyglm(as.formula(pred_demo_formula),
                     design = design_dep, family = quasibinomial())
pred_demo_ors <- extract_or(pred_demo) |>
  mutate(p_adj = p.adjust(p, method = "BH"))

# Model 3: Fully adjusted
pred_formula <- paste("imd_binary ~", paste(model_covars, collapse = " + "))
pred_model <- svyglm(as.formula(pred_formula),
                      design = design_dep, family = quasibinomial())

pred_vif <- tryCatch(car::vif(pred_model), error = function(e) {
  message("  VIF computation failed: ", e$message); NULL
})
if (!is.null(pred_vif)) {
  vif_adj <- if (is.matrix(pred_vif)) pred_vif[, "GVIF^(1/(2*Df))"] else pred_vif
  if (any(vif_adj > sqrt(5)))
    warning("Adjusted VIF > √5 detected: ",
            paste(names(vif_adj[vif_adj > sqrt(5)]), collapse = ", "))
}

pred_ors <- extract_or(pred_model) |>
  mutate(p_adj = p.adjust(p, method = "BH"))

# Model 3 additionally adjusted for harmonised leisure-time physical activity
pred_model_pa <- svyglm(reformulate(model_covars_pa, "imd_binary"),
                        design = design_dep, family = quasibinomial())
pred_ors_pa <- extract_or(pred_model_pa)

pred_model_info <- tibble(
  model       = c("Demographics-adjusted", "Fully adjusted"),
  n           = c(nobs(pred_demo), nobs(pred_model)),
  df_residual = c(pred_demo$df.residual, pred_model$df.residual),
  n_coef      = c(length(coef(pred_demo)), length(coef(pred_model)))
)
model3_cycles <- dep_df |>
  filter(if_all(all_of(model_covars), ~ !is.na(.x))) |>
  count(cycle, name = "n_model3")

# Influence diagnostics for Model 3 (dfbetas, Cook's D)
dfb <- dfbetas(pred_model)
dfb_cut <- 2 / sqrt(nrow(dfb))
cooks <- cooks.distance(pred_model)
influence <- list(
  dfbetas = tibble(
    term = colnames(dfb),
    n_above_cutoff = colSums(abs(dfb) > dfb_cut),
    max_abs_dfbetas = apply(abs(dfb), 2, max)
  ) |> filter(term != "(Intercept)"),
  cutoff = dfb_cut,
  cooks_max = max(cooks),
  cooks_n_above_4_over_n = sum(cooks > 4 / length(cooks)),
  weight_range = range(weights(pred_model)),
  n = nrow(dfb)
)
# Refit Model 3 without observations with Cook's D > 4 / n
model3_rows <- which(complete.cases(design_dep$variables[model_covars]))
drop_seqn <- design_dep$variables$SEQN[model3_rows][cooks > 4 / length(cooks)]
pred_model_trim <- svyglm(reformulate(model_covars, "imd_binary"),
                          design = subset(design_dep,
                                          !SEQN %in% drop_seqn),
                          family = quasibinomial())
influence_refit <- extract_or(pred_model) |>
  select(term, OR, CI_lo, CI_hi) |>
  inner_join(extract_or(pred_model_trim) |>
               select(term, OR_trim = OR, CI_lo_trim = CI_lo,
                      CI_hi_trim = CI_hi),
             by = "term") |>
  mutate(sign_changed = as.integer(sign(log(OR)) != sign(log(OR_trim))))
influence$n_trim <- nobs(pred_model_trim)
influence$n_dropped <- length(drop_seqn)
influence$n_dropped_smallest_groups <- sum(
  design_dep$variables$race_ethnicity[model3_rows][cooks > 4 / length(cooks)] %in%
    c("NH Asian", "Other/Multi-Racial"))
influence$dropped_weight_median_rel <-
  median(weights(pred_model)[cooks > 4 / length(cooks)]) / mean(weights(pred_model))

write_csv(influence$dfbetas |>
            left_join(influence_refit, by = "term") |>
            mutate(cutoff = dfb_cut, cooks_max = influence$cooks_max,
                   n = influence$n,
                   weight_min = influence$weight_range[1],
                   weight_max = influence$weight_range[2],
                   n_dropped = influence$n_dropped,
                   n_trim = influence$n_trim,
                   n_dropped_smallest_groups = influence$n_dropped_smallest_groups,
                   dropped_weight_median_rel = round(influence$dropped_weight_median_rel, 2),
                   across(c(max_abs_dfbetas, cutoff, cooks_max, weight_min,
                            weight_max), ~ round(.x, 3)),
                   across(c(OR, CI_lo, CI_hi, OR_trim, CI_lo_trim, CI_hi_trim),
                          ~ round(.x, 4))),
          here::here("output", "tables", "tableS7_influence_model3.csv"))

# Linearity check for continuous predictors
lin_age <- tryCatch({
  m2 <- svyglm(reformulate(c(model_covars, "I(RIDAGEYR^2)"), "imd_binary"),
                design = design_dep, family = quasibinomial())
  regTermTest(m2, "I(RIDAGEYR^2)")$p[1]
}, error = function(e) NA_real_)

lin_pir <- tryCatch({
  m2 <- svyglm(reformulate(c(model_covars, "I(INDFMPIR^2)"), "imd_binary"),
                design = design_dep, family = quasibinomial())
  regTermTest(m2, "I(INDFMPIR^2)")$p[1]
}, error = function(e) NA_real_)

message("  Linearity: age^2 P = ", round(lin_age, 3),
        ", PIR^2 P = ", round(lin_pir, 3))
message("  Design df (depressed): ", survey::degf(design_dep),
        "; residual df Model 3: ", pred_model$df.residual,
        " (n = ", nobs(pred_model), ")")

# Model diagnostics for the supplement captions
write_csv(tibble(
  age_sq_p = round(lin_age, 3), pir_sq_p = round(lin_pir, 3),
  max_scaled_gvif = round(max(if (is.matrix(pred_vif)) pred_vif[, 3]
                              else pred_vif), 2),
  n_model3 = nobs(pred_model), df_model3 = pred_model$df.residual,
  n_model2 = nobs(pred_demo), df_model2 = pred_demo$df.residual),
  here::here("output", "tables", "model_diagnostics.csv"))

# E-values for significant associations (VanderWeele & Ding 2017)
compute_evalue <- function(or, ci_lo, ci_hi) {
  or <- as.numeric(or); ci_lo <- as.numeric(ci_lo); ci_hi <- as.numeric(ci_hi)
  if (is.na(or)) return(c(e_point = NA_real_, e_ci = NA_real_))
  or_use <- if (or < 1) 1 / or else or
  ci_use <- if (or < 1) 1 / ci_hi else ci_lo
  e_point <- or_use + sqrt(or_use * (or_use - 1))
  e_ci <- if (!is.na(ci_use) && ci_use > 1) {
    ci_use + sqrt(ci_use * (ci_use - 1))
  } else 1
  c(e_point = round(e_point, 2), e_ci = round(e_ci, 2))
}

sig_ors <- pred_ors |> filter(p < 0.05)
evalues <- map_dfr(seq_len(nrow(sig_ors)), \(i) {
  ev <- compute_evalue(sig_ors$OR[i], sig_ors$CI_lo[i], sig_ors$CI_hi[i])
  tibble(term = sig_ors$term[i], e_point = ev["e_point"], e_ci = ev["e_ci"])
})
if (nrow(evalues) > 0)
  message("  E-values: ", paste(evalues$term, "=", evalues$e_point,
                                 "(CI:", evalues$e_ci, ")", collapse = "; "))

# ============================================================
# 5. Race x dep_subtype interaction
# ============================================================
message("  [5/7] Race x subtype interaction...")

interaction_p <- map_dbl(metabolic_vars, \(v) {
  tryCatch({
    m <- svyglm(reformulate(c("dep_subtype", "dep_subtype:race_ethnicity",
                               model_covars),
                             response = v),
                design = design_for(v))
    regTermTest(m, "dep_subtype:race_ethnicity")$p[1]
  }, error = function(e) NA_real_)
}) |> set_names(metabolic_vars)
interaction_p_adj <- bh_adjust(interaction_p)

# ============================================================
# 6. CBC inflammatory validation
# ============================================================
message("  [6/7] CBC inflammatory indices validation...")

cbc_vars <- intersect(c("nlr", "plr", "sii", "wbc"), names(df))
cbc_means <- map(cbc_vars, \(v) {
  tryCatch(
    svyby(reformulate(v), ~dep_subtype, design,
          svymean, na.rm = TRUE),
    error = function(e) NULL
  )
}) |> set_names(cbc_vars) |> compact()

cbc_tests <- map(cbc_vars, \(v) {
  tryCatch(
    svyglm(reformulate("dep_subtype", response = v), design = design),
    error = function(e) NULL
  )
}) |> set_names(cbc_vars) |> compact()

# --- Save results ---
message("  [7/7] Saving results...")
results <- list(
  prevalence = list(
    depression = list(overall = prev_overall, by_race = prev_by_race,
                      by_sex = prev_by_sex, by_age = prev_by_age),
    imd = list(overall = imd_overall, by_race = imd_by_race,
               by_sex = imd_by_sex, by_age = imd_by_age),
    design_df = survey::degf(design)
  ),
  metabolic = list(
    means = metabolic_means, tests = metabolic_tests,
    info = metabolic_info, trig_exam_weight = trig_exam,
    design_df = map_dbl(set_names(metabolic_vars),
                        \(v) survey::degf(design_for(v))),
    p_raw = all_p, p_adj = all_p_adj,
    p_pairwise = pairwise_p, p_pairwise_adj = pairwise_p_adj
  ),
  predictors = list(
    crude = pred_crude, demo = pred_demo_ors,
    model = pred_model, ors = pred_ors, ors_pa = pred_ors_pa, vif = pred_vif,
    model_info = pred_model_info, model3_cycles = model3_cycles,
    influence = influence,
    linearity = list(age_sq_p = lin_age, pir_sq_p = lin_pir),
    evalues = evalues
  ),
  interaction = list(p_raw = interaction_p, p_adj = interaction_p_adj),
  cbc = list(means = cbc_means, tests = cbc_tests),
  cell_sizes = imd_n,
  subgroup_precision = subgroup_precision,
  criterion_overlap = list(table = criterion_overlap, chisq = crit_chisq,
                           diff_model = m_diff, by_definition = symptom_assoc),
  four_group = four_group
)
saveRDS(results, here::here("data", "processed", "analysis_results.rds"))
message("  Results saved to data/processed/analysis_results.rds")
}

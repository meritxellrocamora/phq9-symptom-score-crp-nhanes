source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in  <- here::here("data", "processed",
                   c("analysis_results.rds", "sensitivity_results.rds",
                     "ipcw_mi_results.rds"))
.self <- here::here("R", "04_tables_figures.R")
.out <- here::here("output", "figures", "figS2_flowchart.png")
if (file.exists(.out) &&
    file.mtime(.out) > max(file.mtime(.in)) &&
    file.mtime(.out) > file.mtime(.self)) {
  message("  Up to date, skipping")
} else {

stopifnot(all(file.exists(.in)))
df <- readRDS(here::here("data", "processed", "analytic_dataset.rds"))
results <- readRDS(here::here("data", "processed", "analysis_results.rds"))
sens_results <- readRDS(here::here("data", "processed",
                                   "sensitivity_results.rds"))
ipcw_mi <- readRDS(here::here("data", "processed", "ipcw_mi_results.rds"))
stopifnot(nrow(df) > 0, length(results) > 0)

dir.create(here::here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# Survey design (built once in 02_clean.R); groups carry display labels
design <- readRDS(here::here("data", "processed", "survey_design.rds")) |>
  update(dep_subtype = relabel_groups(dep_subtype))
df <- df |> mutate(dep_subtype = relabel_groups(dep_subtype))

# ============================================================
# Table 1: Demographics by depression subtype
# ============================================================
message("  Table 1...")

tbl1_vars <- c("RIDAGEYR", "sex", "race_ethnicity", "education", "INDFMPIR",
               "bmi", "waist", "smoking", "alcohol", "pa_meets_guidelines",
               "crp", "hba1c", "nlr", "plr", "sii",
               "total_chol", "hdl_chol", "triglycerides")
tbl1_vars <- intersect(tbl1_vars, names(df))

tbl1_labels <- list(
  RIDAGEYR ~ "Age (years)",
  sex ~ "Sex",
  race_ethnicity ~ "Race/Ethnicity",
  education ~ "Education",
  INDFMPIR ~ "Family PIR",
  bmi ~ "BMI (kg/m\u00B2)",
  waist ~ "Waist circumference (cm)",
  smoking ~ "Smoking status",
  alcohol ~ "Alcohol, past 12 months",
  pa_meets_guidelines ~ "Leisure-time MVPA >= 150 min/week",
  crp ~ "CRP (mg/L)",
  hba1c ~ "HbA1c (%)",
  nlr ~ "NLR",
  plr ~ "PLR",
  sii ~ "SII",
  total_chol ~ "Total cholesterol (mg/dL)",
  hdl_chol ~ "HDL cholesterol (mg/dL)"
)
tbl1_digits <- list(all_continuous() ~ 1, all_categorical() ~ 1)

# Later selectors override earlier ones: the median selector comes last
tbl1_main <- design |>
  tbl_svysummary(
    by = dep_subtype,
    include = all_of(setdiff(tbl1_vars, "triglycerides")),
    label = tbl1_labels,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      c(crp, nlr, plr, sii) ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{p}%"
    ),
    digits = tbl1_digits,
    missing = "no"
  ) |>
  add_p() |>
  add_overall()

# Triglycerides: fasting-subsample design (WTSAF2YR / 3, WTSAF2YR > 0)
design_fast <- fasting_design(df)
tbl1_tg <- design_fast |>
  tbl_svysummary(
    by = dep_subtype,
    include = triglycerides,
    label = triglycerides ~ "Triglycerides (mg/dL), fasting subsample\u2020",
    statistic = triglycerides ~ "{median} ({p25}, {p75})",
    digits = tbl1_digits,
    missing = "no"
  ) |>
  add_p() |>
  add_overall()

tbl1 <- tbl_stack(list(tbl1_main, tbl1_tg), quiet = TRUE) |>
  modify_header(all_stat_cols() ~ "**{level}**\nN = {n_unweighted}") |>
  modify_footnote(
    all_stat_cols() ~ "Continuous variables: mean (SD) unless otherwise noted; categorical: %; skewed biomarkers (CRP, NLR, PLR, SII, triglycerides): median (IQR)"
  )

tbl1_gt <- tbl1 |>
  as_gt() |>
  gt::tab_footnote(paste0("Skewed inflammatory and lipid biomarkers (CRP, NLR, ",
                           "PLR, SII, triglycerides) reported as median (IQR); ",
                           "all other continuous variables as mean (SD). ",
                           "Missing data counts (unweighted) are reported ",
                           "in the manuscript text.")) |>
  gt::tab_footnote(paste0("\u2020 Triglycerides are analysed in the fasting ",
                           "subsample with fasting sample weights ",
                           "(WTSAF2YR / 3; participants with WTSAF2YR > 0); ",
                           "the N in the column headers refers to all other ",
                           "rows."))
tbl1_gt |> gt::gtsave(here::here("output", "tables", "table1.html"))
tbl1_gt |> gt::gtsave(here::here("output", "tables", "table1.docx"))

tbl1_csv <- tibble::as_tibble(tbl1)
names(tbl1_csv) <- gsub("\n", " ", gsub("\\*\\*", "", names(tbl1_csv)))
write_csv(tbl1_csv,
          here::here("output", "tables", "table1_characteristics.csv"))
message("  table1.html, table1.docx and table1_characteristics.csv saved")

# ============================================================
# Prevalence estimates table
# ============================================================
message("  Prevalence estimates...")

prev_df <- results$prevalence$design_df
format_prev <- function(svyby_obj, group_var) {
  ci <- confint(svyby_obj, df = prev_df)
  tibble(
    group = rownames(ci),
    prevalence = round(100 * coef(svyby_obj), 1),
    ci_lo = round(100 * ci[, 1], 1),
    ci_hi = round(100 * ci[, 2], 1)
  )
}

prev_table <- bind_rows(
  format_prev(results$prevalence$imd$by_race, "race") |>
    mutate(stratifier = "Race/Ethnicity"),
  format_prev(results$prevalence$imd$by_sex, "sex") |>
    mutate(stratifier = "Sex"),
  format_prev(results$prevalence$imd$by_age, "age") |>
    mutate(stratifier = "Age group")
)
write_csv(prev_table, here::here("output", "tables", "prevalence_estimates.csv"))

# ============================================================
# Metabolic comparison table
# ============================================================
message("  Metabolic comparison table...")

metab_summary <- map2_dfr(results$metabolic$means, names(results$metabolic$means),
  \(svyby_obj, var_name) {
    if (is.null(svyby_obj)) return(NULL)
    tibble(
      marker            = var_name,
      subtype           = as.character(relabel_groups(svyby_obj$dep_subtype)),
      mean              = as.numeric(coef(svyby_obj)),
      se                = as.numeric(SE(svyby_obj)),
      n                 = results$metabolic$info$n[
        match(var_name, results$metabolic$info$marker)],
      df_residual       = results$metabolic$info$df_residual[
        match(var_name, results$metabolic$info$marker)],
      design_df         = results$metabolic$design_df[var_name],
      weights           = if (var_name == "triglycerides")
        "Fasting subsample (WTSAF2YR / 3)"
      else "Examination (2015-2018) / phlebotomy (2021-2023) / 3",
      p_omnibus_raw     = results$metabolic$p_raw[var_name],
      p_omnibus_adj     = results$metabolic$p_adj[var_name],
      p_pairwise_raw    = results$metabolic$p_pairwise[var_name],
      p_pairwise_adj    = results$metabolic$p_pairwise_adj[var_name]
    )
  }
)
fmt_p <- \(p) ifelse(is.na(p), NA_character_,
                     ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
marker_display_names <- c(
  bmi = "BMI (kg/m²)", waist = "Waist circumference (cm)",
  hba1c = "HbA1c (%)", total_chol = "Total cholesterol (mg/dL)",
  hdl_chol = "HDL cholesterol (mg/dL)", triglycerides = "Triglycerides (mg/dL)",
  crp = "CRP (mg/L)", nlr = "NLR", plr = "PLR", sii = "SII"
)
metab_export <- metab_summary |>
  select(-design_df) |>
  mutate(
    mean = round(mean, 2), se = round(se, 3),
    across(starts_with("p_"), fmt_p),
    marker = coalesce(marker_display_names[marker], marker)
  ) |>
  rename(Marker = marker, Subtype = subtype,
         `Survey-weighted mean` = mean, SE = se,
         `Residual df` = df_residual, Weights = weights,
         `p omnibus (unadj)` = p_omnibus_raw, `p omnibus (FDR adj)` = p_omnibus_adj,
         `p pairwise (unadj)` = p_pairwise_raw, `p pairwise (FDR adj)` = p_pairwise_adj)
write_csv(metab_export, here::here("output", "tables", "metabolic_comparison.csv"))

trig_cmp <- bind_rows(
  tibble(weights = "Fasting subsample (WTSAF2YR / 3)",
         subtype = as.character(results$metabolic$means$triglycerides$dep_subtype),
         mean = round(as.numeric(coef(results$metabolic$means$triglycerides)), 2),
         se = round(as.numeric(SE(results$metabolic$means$triglycerides)), 3),
         p_omnibus = results$metabolic$p_raw["triglycerides"],
         p_pairwise_imd_vs_non_imd = results$metabolic$p_pairwise["triglycerides"]),
  tibble(weights = "Examination weight / 3 (previous analysis)",
         subtype = as.character(results$metabolic$trig_exam_weight$means$dep_subtype),
         mean = round(as.numeric(coef(results$metabolic$trig_exam_weight$means)), 2),
         se = round(as.numeric(SE(results$metabolic$trig_exam_weight$means)), 3),
         p_omnibus = results$metabolic$trig_exam_weight$p_omnibus,
         p_pairwise_imd_vs_non_imd = results$metabolic$trig_exam_weight$p_pairwise)
)
write_csv(trig_cmp, here::here("output", "tables",
                               "triglycerides_weight_comparison.csv"))

# ============================================================
# Figure 2: Biomarker profiles (composite A + B; fig2_biomarker_profiles)
#   Panel A: Metabolic markers (dot-whisker)
#   Panel B: CBC inflammatory indices (violin + box)
# ============================================================
message("  Figure 2: biomarker profiles (composite)...")

marker_labels <- c(
  bmi = "BMI (kg/m²)", waist = "Waist (cm)", hba1c = "HbA1c (%)",
  total_chol = "Total cholesterol\n(mg/dL)", hdl_chol = "HDL (mg/dL)",
  triglycerides = "TG (mg/dL)", crp = "CRP (mg/L)",
  nlr = "NLR", plr = "PLR", sii = "SII"
)

marker_labels["crp"] <- "CRP, definitional"

fig2_data <- metab_summary |>
  mutate(
    marker_label = factor(marker_labels[marker], levels = marker_labels),
    subtype = factor(subtype, levels = unname(group_labels)),
    t_crit = qt(0.975, design_df)
  )

fig2a_base <- \(d) {
  ggplot(d, aes(x = subtype, y = mean, color = subtype)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = mean - t_crit * se, ymax = mean + t_crit * se),
                  width = 0.15, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = imd_palette) +
    labs(x = NULL, y = "Survey-weighted mean (95% CI)", color = NULL) +
    theme_minimal(base_family = "sans", base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(face = "bold", size = 8),
          panel.grid.minor = element_blank())
}

# CRP defines IMD (CRP > 1 mg/L): own facet strip, shaded grey
fig2a <- fig2a_base(filter(fig2_data, marker != "crp")) +
  facet_wrap(~marker_label, scales = "free", ncol = 3)
fig2a_crp <- fig2a_base(filter(fig2_data, marker == "crp")) +
  facet_wrap(~marker_label) +
  theme(strip.background = element_rect(fill = "grey80", colour = NA),
        axis.title.y = element_blank())
fig2a <- wrap_elements(full = fig2a + fig2a_crp + plot_layout(widths = c(3, 1)))

cbc_plot_vars <- intersect(c("nlr", "plr", "sii"), names(df))
fig2b_data <- df |>
  select(dep_subtype, all_of(cbc_plot_vars)) |>
  pivot_longer(all_of(cbc_plot_vars), names_to = "marker", values_to = "value") |>
  filter(!is.na(value), !is.na(dep_subtype)) |>
  mutate(marker = toupper(marker))

fig2b <- ggplot(fig2b_data, aes(x = dep_subtype, y = value, fill = dep_subtype)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = imd_palette) +
  labs(x = NULL, y = "Index value (unweighted individual observations)",
       fill = NULL) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank())

fig2_composite <- fig2a / fig2b +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(here::here("output", "figures", "fig2_biomarker_profiles.png"),
       fig2_composite, width = 6.7, height = 8, dpi = 300)
ggsave(here::here("output", "figures", "fig2_biomarker_profiles.pdf"),
       fig2_composite, width = 6.7, height = 8)

# ============================================================
# Figure 1: Sociodemographic patterns (composite A + B; fig1_sociodemographic)
#   Panel A: IMD prevalence by race/ethnicity
#   Panel B: Forest plot, predictors of IMD
# ============================================================
message("  Figure 1: sociodemographic patterns (composite)...")

# Panel A uses the subgroup precision table (logit CIs, unweighted n,
# NCHS presentation standard); estimates failing the standard are hollow
prec <- results$subgroup_precision |> filter(stratifier != "Overall")
fig3a_levels <- c(
  prec |> filter(stratifier == "Race/Ethnicity") |>
    arrange(pct_weighted) |> pull(group),
  prec |> filter(stratifier != "Race/Ethnicity") |> pull(group)
)
fig3a_data <- prec |>
  left_join(prev_table |> select(group, stratifier, ci_lo, ci_hi),
            by = c("group", "stratifier")) |>
  mutate(group = factor(group, levels = fig3a_levels),
         stratifier = factor(stratifier,
                             levels = c("Race/Ethnicity", "Sex", "Age group")),
         meets_nchs = nchs_status == "Meets standard",
         group_n = paste0(group |> str_replace("^NH ", "Non-Hispanic ") |>
                            str_replace("Other/Multi-Racial", "Other or multiracial"),
                          " (", n_imd, "/", n_depressed, ")"),
         group_n = factor(group_n, levels = group_n[order(group)]))

fig3a <- ggplot(fig3a_data, aes(x = pct_weighted, y = group_n)) +
  geom_col(aes(fill = meets_nchs, colour = meets_nchs), width = 0.6) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2,
                orientation = "y") +
  scale_fill_manual(values = c(`TRUE` = "#999999", `FALSE` = "white"),
                    guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "#999999", `FALSE` = "#999999"),
                      guide = "none") +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  facet_grid(stratifier ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Depressed adults meeting the criteria (%)", y = NULL) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(face = "bold", size = 8, angle = 0))

term_labels <- c(
  RIDAGEYR = "Age (per year)", sexFemale = "Female (vs Male)",
  race_ethnicityNH.Black = "Non-Hispanic Black", `race_ethnicityNH Black` = "Non-Hispanic Black",
  race_ethnicityNH.Asian = "Non-Hispanic Asian", `race_ethnicityNH Asian` = "Non-Hispanic Asian",
  `race_ethnicityMexican American` = "Mexican American",
  race_ethnicityMexican.American = "Mexican American",
  `race_ethnicityOther Hispanic` = "Other Hispanic",
  race_ethnicityOther.Hispanic = "Other Hispanic",
  `race_ethnicityOther/Multi-Racial` = "Other or multiracial",
  `race_ethnicityOther.Multi.Racial` = "Other or multiracial",
  educationHS.GED = "High school or equivalent", `educationHS/GED` = "High school or equivalent",
  `educationSome college` = "Some college",
  educationSome.college = "Some college",
  `educationCollege+` = "College or above", educationCollege. = "College or above",
  INDFMPIR = "Income-to-poverty ratio",
  smokingFormer = "Former smoker", smokingCurrent = "Current smoker",
  alcohol = "Alcohol, past 12 months", pa_meets_guidelines = "Leisure-time physical activity"
)

term_groups <- c(
  "Age (per year)" = "Demographics", "Female (vs Male)" = "Demographics",
  "Non-Hispanic Black" = "Race/Ethnicity", "Non-Hispanic Asian" = "Race/Ethnicity",
  "Mexican American" = "Race/Ethnicity", "Other Hispanic" = "Race/Ethnicity",
  "Other or multiracial" = "Race/Ethnicity",
  "High school or equivalent" = "Education", "Some college" = "Education",
  "College or above" = "Education",
  "Income-to-poverty ratio" = "Socioeconomic",
  "Former smoker" = "Lifestyle", "Current smoker" = "Lifestyle",
  "Alcohol, past 12 months" = "Lifestyle", "Leisure-time physical activity" = "Lifestyle"
)

# ============================================================
# Table S6: crude / demographics-adjusted / fully adjusted ORs
# ============================================================
fmt_or <- \(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
s2_part <- function(d, tag) {
  out <- tibble(term = d$term, or_ci = fmt_or(d$OR, d$CI_lo, d$CI_hi),
                p = fmt_p(d$p))
  if ("p_adj" %in% names(d)) out$p_bh <- fmt_p(d$p_adj)
  out <- bind_cols(out, tibble(n = d$n, df = d$df, cycles = d$cycles))
  out |> rename_with(\(x) paste0(tag, "_", x), -term)
}
tableS2 <- s2_part(results$predictors$ors, "adjusted") |>
  left_join(s2_part(results$predictors$crude, "crude"), by = "term") |>
  left_join(s2_part(results$predictors$demo, "demographics"), by = "term") |>
  mutate(predictor = coalesce(term_labels[term], term), .before = 1) |>
  select(predictor, starts_with("crude"), starts_with("demographics"),
         starts_with("adjusted"))
write_csv(tableS2, here::here("output", "tables", "tableS6_sociodemographic.csv"))

# Table S2: included vs excluded (from 02)
incl_excl <- read_csv(here::here("output", "tables", "included_vs_excluded.csv"),
                      show_col_types = FALSE)
write_csv(incl_excl |> rename(Characteristic = variable, Included = included,
                              Excluded = excluded),
          here::here("output", "tables", "tableS2_included_vs_excluded.csv"))

# Table S4: pairwise IMD vs non-IMD depression, non-definitional BH family
s4_wide <- metab_summary |>
  filter(subtype != "No depression") |>
  mutate(cell = sprintf("%.2f (%.2f)", mean, se),
         subtype = if_else(subtype == "Depression, criteria met", "met",
                           "not_met")) |>
  select(marker, subtype, cell, mean, n, p_pairwise_raw, p_pairwise_adj,
         weights) |>
  pivot_wider(names_from = subtype, values_from = c(cell, mean)) |>
  mutate(difference = round(mean_met - mean_not_met, 2),
         Marker = coalesce(marker_display_names[marker], marker),
         definitional = if_else(marker == "crp", "yes", "no")) |>
  transmute(Marker, `Criteria met, mean (SE)` = cell_met,
            `Criteria not met, mean (SE)` = cell_not_met,
            Difference = difference, `n (unweighted, all subtypes)` = n,
            `P (unadjusted)` = fmt_p(p_pairwise_raw),
            `P (BH, 9 non-definitional markers)` = fmt_p(p_pairwise_adj),
            Definitional = definitional, Weights = weights)
write_csv(s4_wide, here::here("output", "tables", "tableS4_pairwise_metabolic.csv"))

# Table S8: complete-case vs IPCW + MI coefficients
s5 <- ipcw_mi$comparison$predictors |>
  mutate(Term = coalesce(term_labels[term], term),
         `Complete-case OR (95% CI)` = fmt_or(OR_cc, CI_lo_cc, CI_hi_cc),
         `Complete-case P` = fmt_p(p_cc),
         `IPCW + MI OR (95% CI)` = fmt_or(OR_mi, CI_lo_mi, CI_hi_mi),
         `IPCW + MI P` = fmt_p(p_mi),
         FMI = round(fmi_coef, 3), `n (CC)` = n_cc, `df (CC)` = df_cc,
         `n (MI)` = n_mi, `df (MI)` = round(df_mi, 1)) |>
  select(Term, `Complete-case OR (95% CI)`, `Complete-case P`,
         `IPCW + MI OR (95% CI)`, `IPCW + MI P`, FMI, `n (CC)`, `df (CC)`,
         `n (MI)`, `df (MI)`)
s5_prev <- ipcw_mi$comparison$prevalence |>
  mutate(Term = paste0("IMD prevalence among depressed: ", analysis),
         `Complete-case OR (95% CI)` = sprintf("%.1f%% (%.1f-%.1f)", imd_prev,
                                              ci_lo, ci_hi)) |>
  select(Term, `Complete-case OR (95% CI)`)
write_csv(bind_rows(s5, s5_prev),
          here::here("output", "tables", "tableS8_ipcw_mi.csv"))

# Software versions for the supplement
pkgs <- c("survey", "mice", "mitools", "gtsummary", "car", "ggplot2",
          "patchwork", "dagitty", "ggdag", "nhanesA", "furrr", "renv")
write_csv(tibble(Software = c("R", pkgs),
                 Version = c(paste(R.version$major, R.version$minor, sep = "."),
                             map_chr(pkgs, \(x) as.character(packageVersion(x))))),
          here::here("output", "tables", "software_versions.csv"))

fig3b_data <- results$predictors$ors |>
  mutate(
    label = coalesce(term_labels[term], term),
    group = factor(term_groups[label],
                   levels = c("Demographics", "Race/Ethnicity", "Education",
                              "Socioeconomic", "Lifestyle")),
    label = fct_reorder2(label, group, OR,
                         .fun = \(g, o) -as.numeric(g[1]) * 100 + o[1])
  )

fig3b <- ggplot(fig3b_data, aes(x = OR, y = label, color = group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), width = 0.25) +
  scale_x_log10() +
  scale_color_manual(values = group_palette) +
  labs(x = "Odds ratio (95% CI)", y = NULL, color = NULL) +
  theme_minimal(base_family = "sans", base_size = 10) +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.width = grid::unit(8, "pt"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 10, 5, 5))

fig3_composite <- (fig3a + theme(plot.margin = margin(5, 15, 5, 5))) + fig3b +
  plot_layout(widths = c(1, 1.6), guides = "collect") &
  theme(legend.position = "bottom")
fig3_composite <- fig3_composite +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(here::here("output", "figures", "fig1_sociodemographic.png"),
       fig3_composite, width = 6.7, height = 5.5, dpi = 300)
ggsave(here::here("output", "figures", "fig1_sociodemographic.pdf"),
       fig3_composite, width = 6.7, height = 5.5)

# ============================================================
# Table S9: sensitivity analyses (primary, 05 analyses, IPCW + MI)
# ============================================================
message("  Table S9: sensitivity analyses...")

pir_summary <- \(or, lo, hi, p, n, d_f) tibble(
  pir_or_ci = if (is.na(or)) NA_character_ else fmt_or(or, lo, hi),
  pir_p = if (is.na(p)) NA_character_ else fmt_p(p),
  model_n = n, model_df = d_f)
prev_summary <- \(prev, d_f = prev_df) tibble(
  imd_prev = round(100 * as.numeric(coef(prev)), 1),
  imd_ci = sprintf("(%.1f-%.1f)", 100 * confint(prev, df = d_f)[1],
                   100 * confint(prev, df = d_f)[2]))

pir_row <- results$predictors$ors |> filter(term == "INDFMPIR")
s6_primary <- bind_cols(
  tibble(analysis = "Primary (complete case)", n_total = nrow(df),
         n_depressed = sum(df$depression == 1),
         n_imd = sum(df$imd_binary == 1)),
  prev_summary(results$prevalence$imd$overall),
  pir_summary(pir_row$OR, pir_row$CI_lo, pir_row$CI_hi, pir_row$p,
              pir_row$n, pir_row$df))
s6_sens <- map_dfr(sens_results, \(x) {
  pm <- x$predictors
  pir <- if (is.null(pm)) NULL else {
    cf <- summary(pm)$coefficients["INDFMPIR", ]
    crit <- qt(0.975, pm$df.residual)
    list(or = exp(cf[["Estimate"]]),
         lo = exp(cf[["Estimate"]] - crit * cf[["Std. Error"]]),
         hi = exp(cf[["Estimate"]] + crit * cf[["Std. Error"]]),
         p = cf[["Pr(>|t|)"]], n = nobs(pm), df = pm$df.residual)
  }
  bind_cols(
    tibble(analysis = x$label, n_total = x$n_total, n_depressed = x$n_dep,
           n_imd = x$n_imd),
    prev_summary(x$prevalence, x$design_df),
    if (is.null(pir)) pir_summary(NA, NA, NA, NA, NA_integer_, NA_real_)
    else pir_summary(pir$or, pir$lo, pir$hi, pir$p, pir$n, pir$df))
})
ipcw_pir <- ipcw_mi$pooled$predictors |> filter(term == "INDFMPIR")
ipcw_prev <- ipcw_mi$pooled$prev_imd
s6_ipcw <- bind_cols(
  tibble(analysis = "IPCW + MI (m = 20)", n_total = ipcw_mi$n$analytic,
         n_depressed = sum(df$depression == 1), n_imd = sum(df$imd_binary == 1),
         imd_prev = round(100 * ipcw_prev$estimate, 1),
         imd_ci = sprintf("(%.1f-%.1f)", 100 * ipcw_prev$CI_lo,
                          100 * ipcw_prev$CI_hi)),
  pir_summary(ipcw_pir$OR, ipcw_pir$CI_lo, ipcw_pir$CI_hi, ipcw_pir$p,
              ipcw_pir$n_model, round(ipcw_pir$df, 1)))
tableS6 <- bind_rows(s6_primary, s6_sens, s6_ipcw)
write_csv(tableS6, here::here("output", "tables", "tableS9_sensitivity.csv"))

# ============================================================
# Flowchart (grid graphics)
# ============================================================
message("  Flowchart...")

flow <- read_csv(here::here("output", "tables", "exclusion_flowchart.csv"),
                 show_col_types = FALSE)

draw_flowchart <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = 0.9, height = 0.9))

  box <- function(label, x, y, w = 0.4, h = 0.08) {
    grid::grid.rect(x = x, y = y, width = w, height = h,
                    gp = grid::gpar(fill = "white", col = "black"))
    grid::grid.text(label, x = x, y = y,
                    gp = grid::gpar(fontsize = 9))
  }
  arrow <- function(x1, y1, x2, y2) {
    grid::grid.lines(x = c(x1, x2), y = c(y1, y2),
                     arrow = grid::arrow(length = grid::unit(0.02, "npc")))
  }

  n <- flow
  fmt <- \(x) format(x, big.mark = ",")
  ys <- c(0.92, 0.72, 0.52, 0.32, 0.12)
  box(paste0("Adults aged 18 years or older interviewed,\n",
             "NHANES 2015-2016, 2017-2018, 2021-2023\nN = ", fmt(n$n_adults)),
      0.4, ys[1], w = 0.5, h = 0.13)
  arrow(0.4, ys[1] - 0.065, 0.4, ys[2] + 0.05)
  box(paste0("Aged 20 years or older, non-pregnant, with examination\n",
             "weight (2015-2018) or phlebotomy weight (2021-2023)\nN = ",
             fmt(n$n_total)),
      0.4, ys[2], w = 0.5, h = 0.1)
  arrow(0.4, ys[2] - 0.05, 0.4, ys[3] + 0.04)
  box(paste0("Complete PHQ-9\nN = ", fmt(n$n_phq9)), 0.4, ys[3], w = 0.5)
  arrow(0.4, ys[3] - 0.04, 0.4, ys[4] + 0.04)
  box(paste0("CRP measured and 10 mg/L or below\nN = ", fmt(n$n_crp)),
      0.4, ys[4], w = 0.5)
  arrow(0.4, ys[4] - 0.04, 0.4, ys[5] + 0.04)
  box(paste0("Complete age, sex, race/ethnicity, education\nand ",
             "income-to-poverty ratio: analytic sample\nN = ",
             fmt(n$n_covars)), 0.4, ys[5], w = 0.5, h = 0.1)

  # Exclusion boxes (right side) with connecting arrows
  ex <- \(label, y, h = 0.08) {
    arrow(0.4, y, 0.66, y)
    box(label, 0.82, y, w = 0.3, h = h)
  }
  ex(paste0("Excluded: aged 18-19 years, n = ", fmt(n$n_age_18_19),
            "\nzero weight, n = ",
            fmt(n$n_adults - n$n_age_18_19 - n$n_weighted),
            "\npregnant, n = ", fmt(n$n_pregnant)),
     mean(ys[1:2]), h = 0.11)
  ex(paste0("Excluded: PHQ-9 incomplete\nn = ", fmt(n$n_total - n$n_phq9)),
     mean(ys[2:3]))
  ex(paste0("Excluded: CRP missing, n = ", fmt(n$n_crp_missing),
            "\nCRP above 10 mg/L (acute), n = ", fmt(n$n_crp_acute)),
     mean(ys[3:4]))
  ex(paste0("Excluded: covariates missing\nn = ", fmt(n$n_crp - n$n_covars)),
     mean(ys[4:5]))
}

tryCatch({
  pdf(here::here("output", "figures", "figS2_flowchart.pdf"),
      width = 8, height = 6)
  draw_flowchart()
}, finally = dev.off())

tryCatch({
  png(here::here("output", "figures", "figS2_flowchart.png"),
      width = 8, height = 6, units = "in", res = 300)
  draw_flowchart()
}, finally = dev.off())

message("  All tables and figures saved.")
}

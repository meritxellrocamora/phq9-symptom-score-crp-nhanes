cycles <- tibble::tibble(
  id       = c("I", "J", "L"),
  label    = c("2015-2016", "2017-2018", "2021-2023"),
  type     = rep("suffix", 3),
  duration = rep(2L, 3),
  # WTPH2YR exists only from 2021-2023; HSCRP_I/J carry no subsample weight
  wt_blood = c("WTMEC2YR", "WTMEC2YR", "WTPH2YR")
)
n_cycles <- nrow(cycles)

# IMD thresholds (Lamers 2024)
crp_cutoff_primary <- 1
crp_cutoff_strict  <- 3
crp_acute_cutoff   <- 10
aes_threshold      <- 6
depression_cutoff  <- 10

# CRP bridging to the 2017-2018 Cobas 6000 scale (NCHS analytic notes):
#   2015-2016 Beckman DxC 660i -> Cobas 6000, forward equation in HSCRP_J,
#   valid for DxC values <= 23 mg/L (higher values left unadjusted);
#   2021-2023 Cobas 8000 -> Cobas 6000, backward equation in HSCRP_L.
crp_xwalk_I_slope <- 0.8695
crp_xwalk_I_int   <- 0.2954
crp_xwalk_I_max   <- 23
crp_xwalk_int   <- 0.08946
crp_xwalk_slope <- 1.034
crp_log_offset  <- 0.01

# PHQ-9 and AES item vectors
phq9_vars <- paste0("DPQ0", 1:9, "0")
aes_items <- c("DPQ030", "DPQ040", "DPQ050", "DPQ080")

# Model covariates (Models 1-3, sensitivity, IPCW + MI). Physical activity
# is excluded from the primary adjustment set: the 2021-2023 questionnaire
# replaced the GPAQ items with leisure-time frequency/duration items without
# a 10-minute bout threshold, and the harmonised leisure-time variable
# shifts from 40-42% (2015-2018) to 57% (2021-2023) in the analytic sample,
# an instrument effect.
# It is kept as a descriptive variable and in sensitivity model S11.
model_covars <- c("RIDAGEYR", "sex", "race_ethnicity", "education",
                  "INDFMPIR", "smoking", "alcohol")
model_covars_pa <- c(model_covars, "pa_meets_guidelines")
# Covariates with a documented measurement discontinuity between cycles
# (drift check reports them with message() rather than warning())
known_cycle_discontinuity <- c("pa_meets_guidelines")

# Fasting-subsample design for triglycerides (WTSAF2YR / n_cycles)
fasting_design <- function(data) {
  svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_fasting,
            nest = TRUE, data = data) |>
    subset(wt_fasting > 0)
}

# Display labels for the three groups (internal levels keep the short names)
group_labels <- c("No depression" = "No depression",
                  "Non-IMD depression" = "Depression, criteria not met",
                  "IMD" = "Depression, criteria met")
relabel_groups <- \(x) factor(group_labels[as.character(x)],
                              levels = unname(group_labels))

# Palette
imd_palette <- c("No depression" = "#88CCEE",
                 "Depression, criteria not met" = "#DDCC77",
                 "Depression, criteria met" = "#CC6677")

# Figure 1B term groups (Okabe-Ito)
group_palette <- c("Demographics" = "#E69F00", "Race/Ethnicity" = "#56B4E9",
                   "Education" = "#009E73", "Socioeconomic" = "#D55E00",
                   "Lifestyle" = "#CC79A7")

# Reproducibility
set.seed(12345)

source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in  <- here::here("R", "02b_dag.R")
.out <- here::here("output", "figures", "figS1_dag.pdf")
if (file.exists(.out) && file.mtime(.out) > file.mtime(.in)) {
  message("  Up to date, skipping")
} else {

message("  Specifying DAG and adjustment sets...")

# Nodes:
#   IMD classification (outcome): defined by CRP + AES among depressed
#   Confounders: age, sex, race/ethnicity, education, PIR, smoking,
#                alcohol, physical activity
#   Mediators: BMI, metabolic markers (on causal path -> do NOT adjust)
#   Mediator: healthcare access (on causal path -> do NOT adjust)

dag <- dagitty::dagitty('
dag {
  IMD [outcome]
  Age -> IMD
  Age -> BMI
  Age -> MetabolicMarkers
  Sex -> IMD
  Sex -> BMI
  Race -> IMD
  Race -> BMI
  Race -> MetabolicMarkers
  Education -> IMD
  Education -> BMI
  Education -> Income
  Education -> HealthcareAccess
  Income -> IMD
  Income -> BMI
  Income -> HealthcareAccess
  Age -> Smoking
  Age -> Alcohol
  Age -> PhysicalActivity
  Sex -> Smoking
  Sex -> Alcohol
  Sex -> PhysicalActivity
  Smoking -> IMD
  Smoking -> BMI
  Smoking -> MetabolicMarkers
  Alcohol -> IMD
  Alcohol -> MetabolicMarkers
  PhysicalActivity -> IMD
  PhysicalActivity -> BMI
  BMI -> MetabolicMarkers
  BMI -> IMD
  MetabolicMarkers -> IMD
  Race -> HealthcareAccess
  Race -> Education
  Race -> Income
  Race -> Smoking
  HealthcareAccess -> IMD
}
')

# Minimally sufficient adjustment set (excluding mediators BMI/Metabolic)
adj_set <- dagitty::adjustmentSets(dag, exposure = "Race", outcome = "IMD",
                                   type = "minimal")
message("Minimally sufficient adjustment sets for Race -> IMD:")
print(adj_set)

# For the predictors model: adjust for confounders only
# Do NOT include BMI or MetabolicMarkers (mediators)
# Do NOT include HealthcareAccess (mediator)
message("\nPredictors model covariates (from DAG):")
message("  Age, Sex, Race, Education, Income, Smoking, Alcohol, PhysicalActivity")

# --- Plot DAG ---
dir.create(here::here("output", "figures"), recursive = TRUE,
           showWarnings = FALSE)

dag_coords <- list(
  x = c(Age = 0, Sex = 1, Race = 3, Education = 0.5, Income = 1.5,
        Smoking = 2.5, Alcohol = 3.5, PhysicalActivity = 4.5,
        BMI = 1.5, MetabolicMarkers = 3, IMD = 2.25,
        HealthcareAccess = 4.5),
  y = c(Age = 2, Sex = 2, Race = 2, Education = 1, Income = 1,
        Smoking = 1, Alcohol = 1, PhysicalActivity = 1,
        BMI = 0, MetabolicMarkers = 0, IMD = -1,
        HealthcareAccess = 0)
)
dagitty::coordinates(dag) <- dag_coords

# Physical activity is measured but kept out of the primary adjustment set
# (instrument change in 2021-2023; sensitivity model S11)
node_types <- c(
  IMD = "Outcome", BMI = "Intermediate (not adjusted)",
  MetabolicMarkers = "Intermediate (not adjusted)",
  HealthcareAccess = "Intermediate (not adjusted)",
  Age = "Adjusted", Sex = "Adjusted", Race = "Adjusted",
  Education = "Adjusted", Income = "Adjusted", Smoking = "Adjusted",
  Alcohol = "Adjusted", PhysicalActivity = "Sensitivity model only"
)

node_labels <- c(
  Age = "Age", Sex = "Sex", Race = "Race", Education = "Education",
  Income = "Income", Smoking = "Smoking", Alcohol = "Alcohol",
  PhysicalActivity = "Physical\nActivity", BMI = "BMI",
  MetabolicMarkers = "Metabolic\nMarkers", IMD = "Criteria\nmet",
  HealthcareAccess = "Healthcare\nAccess"
)

dag_tidy <- tidy_dagitty(dag) |>
  mutate(node_type = node_types[name],
         label = node_labels[name])

dag_plot <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey50") +
  geom_dag_point(aes(color = node_type), size = 15, alpha = 0.9) +
  geom_dag_text(aes(label = label), size = 2.6, color = "black") +
  scale_color_manual(
    values = c("Adjusted" = "#88CCEE", "Sensitivity model only" = "#999933",
               "Intermediate (not adjusted)" = "#DDCC77",
               "Outcome" = "#CC6677"),
    na.value = "grey80"
  ) +
  labs(color = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 5))) +
  theme_dag() +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

ggsave(here::here("output", "figures", "figS1_dag.pdf"),
       dag_plot, width = 6.7, height = 5, bg = "white")
ggsave(here::here("output", "figures", "figS1_dag.png"),
       dag_plot, width = 6.7, height = 5, dpi = 300, bg = "white")

message("  DAG saved: output/figures/figS1_dag.pdf, figS1_dag.png")
}

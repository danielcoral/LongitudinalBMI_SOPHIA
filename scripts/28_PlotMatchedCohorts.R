# Rebuild figures from saved estimates without loading or refitting models.
source("src/project_setup.R")
project_dir <- getwd()
source(file.path(project_dir, "src", "matched_cohort_helpers.R"))
contrasts <- readRDS(file.path(
  project_paths$cross_cohort, "matched_cohort_contrasts.rds"
))
dir.create(file.path(project_paths$plots), showWarnings = FALSE)
for (reference in c("ALSPAC", "RS")) {
  figure <- matched_cohort_plot(contrasts, reference)
  ggplot2::ggsave(
    file.path(project_paths$plots, paste0(
      "matched_", tolower(reference), "_ukb.png"
    )),
    plot = figure, width = 7, height = 4.5,
    units = "in", dpi = 300, bg = "white"
  )
}

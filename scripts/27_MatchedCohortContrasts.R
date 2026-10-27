# Run from the project root after receiving all three cohort result files.
library(mgcv)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "matched_cohort_helpers.R"))
cohorts <- c("ALSPAC", "RS", "UKB")
input_files <- file.path(
  project_dir, "data", cohorts, paste0("resultdf_", cohorts, ".rds")
)
if (!all(file.exists(input_files))) {
  stop("Missing cohort result files: ",
       paste(input_files[!file.exists(input_files)], collapse = ", "))
}
result_table <- purrr::map2_dfr(
  input_files, cohorts,
  function(path, cohort_name) {
    result <- readRDS(path)
    if (!"reftime" %in% names(result)) result$reftime <- NA_character_
    mutate(result, cohort = cohort_name)
  }
)
contrasts <- matched_cohort_contrasts(result_table)
output_dir <- file.path(project_dir, "data", "cross_cohort")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(contrasts, file.path(output_dir, "matched_cohort_contrasts.rds"))
readr::write_csv(
  contrasts, file.path(output_dir, "matched_cohort_contrasts.csv")
)
if (any(!contrasts$within_ukb_central_range)) {
  warning("Some matched values are outside central UKB marginal ranges; ",
          "inspect the support flags before interpreting these estimates.")
}

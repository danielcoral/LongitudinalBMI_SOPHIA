# Within-cohort contrasts: upper versus lower quartile of recent BMI slope at
# the same BMI, in ALSPAC, UK Biobank and the Rotterdam Study.
# Inputs: data/{ALSPAC,RS,UKB}/resultdf_*.rds (saved fitted models).
# Outputs: data/cross_cohort/within_cohort_contrasts.{rds,csv}.
library(mgcv)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "matched_cohort_helpers.R"))
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
results <- load_cohort_results(project_dir)
contrasts <- within_cohort_contrasts(results)

# Checks against earlier outputs; a silent error here would change the paper.
stopifnot(nrow(contrasts) == 318)
earlier <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "matched_cohort_contrasts.csv"),
  show_col_types = FALSE
) |>
  filter(reference_cohort == estimate_cohort) |>
  transmute(
    cohort = reference_cohort, sex, trait, bmi_percentile,
    old_estimate = estimate, old_se = std_error
  )
compared <- inner_join(
  contrasts, earlier, by = c("cohort", "sex", "trait", "bmi_percentile")
)
stopifnot(
  nrow(compared) == 210,
  max(abs(compared$estimate - compared$old_estimate)) < 1e-8,
  max(abs(compared$std_error - compared$old_se)) < 1e-8
)
# The UK Biobank outcome SD must equal 10 x the threshold (0.1 SD) that
# script 34 stored for the +-0.5 kg/m2/year contrasts.
thresholds <- readr::read_csv(
  file.path(project_dir, "data", "UKB", "results",
            "slope_reference_contrasts.csv"),
  show_col_types = FALSE
) |>
  filter(model_family %in% c("primary", "imaging")) |>
  distinct(sex, trait, n_old = n, thresh)
sd_check <- contrasts |>
  filter(cohort == "UKB", bmi_percentile == "p50") |>
  inner_join(thresholds, by = c("sex", "trait"))
stopifnot(
  nrow(sd_check) == 36,
  max(abs(sd_check$outcome_sd / (10 * sd_check$thresh) - 1)) < 1e-10,
  all(sd_check$n == sd_check$n_old)
)

output_dir <- file.path(project_dir, "data", "cross_cohort")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(contrasts, file.path(output_dir, "within_cohort_contrasts.rds"))
readr::write_csv(
  contrasts, file.path(output_dir, "within_cohort_contrasts.csv")
)

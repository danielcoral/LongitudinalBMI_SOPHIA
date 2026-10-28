# Gain, stable and loss contrasts (slopes +0.5, 0 and -0.5 kg/m2/year) and
# the gain-loss asymmetry at one point per cohort: BMI 30 kg/m2 at age 60 in
# UK Biobank (midlife) and age 70 in the Rotterdam Study (older age), close
# to each cohort's median age. Covers the 17 traits measured in both cohorts
# plus log liver fat (UK Biobank only). Other covariates are at each model's
# reference values, as in the within-cohort contrasts (script 36). ALSPAC is
# not included: zero and -0.5 kg/m2/year lie in the extreme lower tail of its
# slope distribution.
# Inputs: data/{UKB,RS}/resultdf_{cohort}.rds (saved models, evaluated here);
#   data/UKB/processed_data/slopebmregdf_2ilf.rds (log liver-fat models);
#   data/cross_cohort/slope_support_summary.csv (script 37).
# Output: stable_contrasts_common_grid.csv in data/cross_cohort/.
library(mgcv)
library(dplyr, warn.conflicts = FALSE)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "matched_cohort_helpers.R"))
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "slope_contrast_helpers.R"))
bmi_value <- 30
ages <- c(UKB = 60, RS = 70)

liver <- readRDS(
  file.path(project_dir, "data", "UKB", "processed_data",
            "slopebmregdf_2ilf.rds")
) |>
  filter(trait == "logLiverF_perc") |>
  transmute(cohort = "UKB", reftime = "2i", sex, trait, mod1)
models <- load_cohort_results(project_dir) |>
  filter(cohort %in% names(ages)) |>
  cohort_models() |>
  select(cohort, reftime, sex, trait, mod1) |>
  bind_rows(liver)

rows <- pmap_dfr(
  models[c("cohort", "reftime", "sex", "trait", "mod1")],
  function(cohort, reftime, sex, trait, mod1) {
    grid <- model_grid(mod1, ages[[cohort]], bmi_value)
    paired_slope_contrasts(mod1, grid) |>
      transmute(
        cohort = cohort, reftime = reftime, sex = sex, trait = trait,
        age, bmi, comparison, estimate, std_error, conf_low, conf_high,
        p_value,
        ratio = if_else(trait == "logLiverF_perc", exp(estimate), NA_real_),
        outcome_sd = outcome_sd(mod1), n = model_n(mod1)
      )
  }
)

# Is the chosen age and BMI inside the central 95% of the cohort's sample?
support <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "slope_support_summary.csv"),
  show_col_types = FALSE
) |>
  select(cohort, reftime, sex, age_p025, age_p975, bmi_p025, bmi_p975)
common <- rows |>
  left_join(support, by = c("cohort", "reftime", "sex")) |>
  left_join(trait_groups, by = "trait") |>
  mutate(
    group = coalesce(group, "Liver fat"),
    estimate_sd = estimate / outcome_sd,
    conf_low_sd = conf_low / outcome_sd,
    conf_high_sd = conf_high / outcome_sd,
    age_in_central95 = age >= age_p025 & age <= age_p975,
    bmi_in_central95 = bmi >= bmi_p025 & bmi <= bmi_p975
  ) |>
  select(
    cohort, reftime, sex, trait, group, age, bmi, comparison, estimate,
    std_error, conf_low, conf_high, p_value, ratio, outcome_sd, estimate_sd,
    conf_low_sd, conf_high_sd, n, age_in_central95, bmi_in_central95
  )
# 17 traits x 2 sexes x 4 comparisons in Rotterdam; 18 traits in UKB.
stopifnot(
  nrow(common) == 136 + 144, !anyNA(common$estimate),
  !anyNA(common$group), all(common$age_in_central95),
  all(common$bmi_in_central95)
)

output_dir <- file.path(project_dir, "data", "cross_cohort")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  common, file.path(output_dir, "stable_contrasts_common_grid.csv")
)

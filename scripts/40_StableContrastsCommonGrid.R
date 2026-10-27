# Gain, stable and loss contrasts (slopes +0.5, 0 and -0.5 kg/m2/year) and
# the gain-loss asymmetry on one grid for UK Biobank and the Rotterdam Study:
# age 65 and BMI 22.5, 27.5 and 32.5 kg/m2, for the 17 traits measured in both
# cohorts plus log liver fat (UK Biobank only). ALSPAC is not included: zero
# kg/m2/year lies in the extreme lower tail of its slope distribution.
# Inputs: data/RS/resultdf_RS.rds (saved models, evaluated here);
#   data/UKB/results/slope_reference_contrasts.csv (script 34, read as is);
#   data/cross_cohort/slope_support_summary.csv (script 38).
# Output: stable_contrasts_common_grid.csv in data/cross_cohort/.
library(mgcv)
library(dplyr, warn.conflicts = FALSE)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "matched_cohort_helpers.R"))
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "slope_contrast_helpers.R"))
bmi_grid <- c(22.5, 27.5, 32.5)

# Rotterdam: the saved models (same contrast function as script 34).
rs_models <- load_cohort_results(project_dir) |>
  filter(cohort == "RS") |>
  cohort_models()
rs_rows <- pmap_dfr(
  rs_models[c("cohort", "reftime", "sex", "trait", "mod1")],
  function(cohort, reftime, sex, trait, mod1) {
    paired_slope_contrasts(mod1, model_grid(mod1, 65, bmi_grid)) |>
      transmute(
        cohort = cohort, reftime = reftime, sex = sex, trait = trait,
        age, bmi, comparison, estimate, std_error, conf_low, conf_high,
        p_value, ratio = NA_real_, outcome_sd = outcome_sd(mod1),
        n = model_n(mod1)
      )
  }
)

# UK Biobank: script 34 already evaluated this grid. Raw liver fat is
# replaced by log liver fat; UK Biobank outcome SD = 10 x the 0.1 SD
# threshold that script 34 stored.
ukb_rows <- readr::read_csv(
  file.path(project_dir, "data", "UKB", "results",
            "slope_reference_contrasts.csv"),
  show_col_types = FALSE
) |>
  filter(age == 65, bmi %in% bmi_grid, trait != "LiverF_perc") |>
  transmute(
    cohort = "UKB", reftime = if_else(model_family == "primary", "rec", "2i"),
    sex, trait, age, bmi, comparison, estimate, std_error, conf_low,
    conf_high, p_value, ratio, outcome_sd = 10 * thresh, n
  )

# Is age 65 and each BMI inside the central 95% of the cohort's sample?
support <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "slope_support_summary.csv"),
  show_col_types = FALSE
) |>
  select(cohort, reftime, sex, age_p025, age_p975, bmi_p025, bmi_p975)
common <- bind_rows(rs_rows, ukb_rows) |>
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
# 17 traits x 2 sexes x 3 BMI x 4 comparisons in Rotterdam; 18 traits in UKB.
stopifnot(
  nrow(common) == 408 + 432, !anyNA(common$estimate),
  !anyNA(common$group), !anyNA(common$age_in_central95)
)

output_dir <- file.path(project_dir, "data", "cross_cohort")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  common, file.path(output_dir, "stable_contrasts_common_grid.csv")
)

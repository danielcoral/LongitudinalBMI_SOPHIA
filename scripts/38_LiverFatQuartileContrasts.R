# Liver fat (log scale): upper versus lower quartile of recent BMI slope.
# The slope pair is the sex-specific 25th and 75th percentile of the slope in
# the liver-MRI sample; ratio = exp(estimate) = upper / lower quartile.
# Nothing is refitted: the saved models are evaluated.
# Input: data/UKB/processed_data/slopebmregdf_2ilf.rds (log liver-fat models
#   with their data).
# Adds age-specific contrasts to the saved liver-fat output.
library(mgcv)
library(dplyr, warn.conflicts = FALSE)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "slope_contrast_helpers.R"))
source(file.path(project_dir, "src", "liver_fat_interaction_helpers.R"))
source(file.path(project_dir, "src", "liver_fat_quartile_helpers.R"))
processed_dir <- file.path(project_dir, "data", "UKB", "processed_data")
results_dir <- file.path(project_dir, "data", "UKB", "results")
liver <- readRDS(file.path(processed_dir, "slopebmregdf_2ilf.rds")) |>
  filter(trait == "logLiverF_perc")
stopifnot(nrow(liver) == 2)
slope_pairs <- map(
  set_names(liver$bmdf, liver$sex),
  function(data) as.numeric(quantile(data$Slope, c(0.25, 0.75)))
)

# "cohort_percentile" evaluates each model at its own sample's BMI
# quartiles and median age (as for every other trait in script 36);
# "ukb_fixed" retains the fixed grid of script 33 for existing summaries.
# "age_bmi_percentile" evaluates MRI-sample BMI quartiles at ages 55, 60 and
# 65 for Figure 3.
liver_rows <- pmap_dfr(
  list(sex = liver$sex, data = liver$bmdf, model = liver$mod1),
  function(sex, data, model) {
    bmi_quartiles <- as.numeric(quantile(data$bmi, c(0.25, 0.5, 0.75)))
    percentile <- quartile_contrast_rows(
      model, model_grid(model, median(data$age), bmi_quartiles),
      slope_pairs[[sex]]
    )
    percentile$grid <- "cohort_percentile"
    percentile$bmi_percentile <- c("p25", "p50", "p75")
    fixed <- quartile_contrast_rows(
      model, slope_contrast_grid(), slope_pairs[[sex]]
    )
    fixed$grid <- "ukb_fixed"
    age_bmi <- model_grid(
      model,
      age = rep(c(55, 60, 65), each = 3),
      bmi = rep(bmi_quartiles, times = 3)
    ) |>
      mutate(
        bmi_percentile = rep(c("p25", "p50", "p75"), times = 3)
      )
    age_bmi <- quartile_contrast_rows(
      model, age_bmi, slope_pairs[[sex]]
    )
    age_bmi$grid <- "age_bmi_percentile"
    age_limits <- quantile(data$age, c(0.025, 0.975))
    bmi_limits <- quantile(data$bmi, c(0.025, 0.975))
    bind_rows(percentile, fixed, age_bmi) |>
      mutate(
        cohort = "UKB", reftime = "2i", sex = sex, trait = "logLiverF_perc",
        p_value = 2 * pnorm(-abs(estimate / std_error)),
        outcome_sd = outcome_sd(model), n = model_n(model),
        estimate_sd = estimate / outcome_sd,
        conf_low_sd = conf_low / outcome_sd,
        conf_high_sd = conf_high / outcome_sd,
        age_in_central95 = age >= age_limits[1] & age <= age_limits[2],
        bmi_in_central95 = bmi >= bmi_limits[1] & bmi <= bmi_limits[2]
      )
  }
) |>
  select(
    cohort, reftime, sex, trait, grid, bmi_percentile, age, bmi, slope_p25,
    slope_p75, estimate, std_error, conf_low, conf_high, p_value, ratio,
    ratio_low, ratio_high, outcome_sd, estimate_sd, conf_low_sd,
    conf_high_sd, n, age_in_central95, bmi_in_central95
  )
stopifnot(nrow(liver_rows) == 48, !anyNA(liver_rows$estimate))

readr::write_csv(
  liver_rows, file.path(results_dir, "liver_fat_quartile_contrasts.csv")
)

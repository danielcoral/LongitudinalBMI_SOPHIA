# Does waist-to-hip ratio (WHR) modify the slope/liver-fat association?
# Same BMI and same recent BMI slope, but central (higher WHR) versus more
# peripheral (lower WHR) fat: does liver fat differ, and does the slope
# contrast differ? WHR is measured at the imaging visit, as is liver fat.
# Inputs: data/UKB/processed_data/slopebmregdf_2ilf.rds (log liver-fat model
#   data); data/UKB/raw_data/WaistHip_Instance2.tsv (UK Biobank fields 48 and
#   49, instance 2; extracted from the main UK Biobank table with
#   `cut -f1,75,79 41806.tab`).
# Outputs (data/UKB/results/): liver_fat_whr_tests.csv,
#   liver_fat_whr_diagnostics.csv, liver_fat_whr_profiles.csv,
#   liver_fat_whr_predictions.csv, liver_fat_whr_slope_contrasts.csv,
#   liver_fat_whr_profile_contrasts.csv,
#   liver_fat_whr_profile_differences.csv.
# Model output: data/UKB/processed_data/liver_fat_whr_models.rds
# dplyr resolves column names through non-standard evaluation.
# nolint start: object_usage_linter.
library(mgcv)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
source("src/liver_fat_interaction_helpers.R")

results_dir <- "data/UKB/results"
model_path <- "data/UKB/processed_data/liver_fat_whr_models.rds"
bmi_values <- c(25, 30, 35)
profile_age <- 60
bmi_half_width <- 2.5
whr_half_width <- 0.025
age_half_width <- 5
profile_levels <- c("Lower (p25)", "Median (p50)", "Higher (p75)")

original <- readRDS("data/UKB/processed_data/slopebmregdf_2ilf.rds") |>
  filter(trait == "logLiverF_perc")
waist_hip <- read_tsv(
  "data/UKB/raw_data/WaistHip_Instance2.tsv", show_col_types = FALSE
) |>
  transmute(eid = f.eid, WHR = f.48.2.0 / f.49.2.0) |>
  filter(is.finite(WHR))
stopifnot(nrow(original) == 2L, !anyDuplicated(waist_hip$eid))

covars <- c(
  "CurrSmoke", "Diabetes", "Insulin", "LipidLower", "HTMed",
  "TSI", "CollegeEd", "Menopause", "ancestry"
)
base_terms <- c(
  covars, "s(age, bs = 'cr')", "s(bmi, bs = 'cr')",
  "ti(age, bmi)", "s(Slope, bs = 'cr')", "ti(bmi, Slope)",
  "ti(Slope, age)", "s(WHR, bs = 'cr')"
)
reduced_formula <- reformulate(base_terms, response = "value")
full_formula <- reformulate(
  c(base_terms, "ti(Slope, WHR)"), response = "value"
)

# Two log-liver contrasts, each between two rows that differ only as stated.
contrast_rows <- function(model, high, low, info) {
  bind_cols(info, log_liver_contrast(model, high, low))
}

fits <- tests <- diagnostics <- profiles <- predictions <- list()
slope_contrasts <- profile_contrasts <- differences <- list()
for (i in seq_len(nrow(original))) {
  sex_value <- original$sex[i]
  # Slope quartiles of the full liver-MRI sample, as in script 39.
  pair <- as.numeric(quantile(original$bmdf[[i]]$Slope, c(0.25, 0.75)))
  slope_median <- median(original$bmdf[[i]]$Slope)
  dat <- inner_join(original$bmdf[[i]], waist_hip, by = "eid")
  stopifnot(!anyDuplicated(dat$eid), all(complete.cases(dat)))
  cat("Fitting", sex_value, "N =", nrow(dat), "\n")
  # fREML matches the original estimation; ML supplies a nested comparison.
  full <- bam(full_formula, data = dat, method = "fREML")
  reduced <- bam(reduced_formula, data = dat, method = "fREML")
  full_ml <- bam(full_formula, data = dat, method = "ML")
  reduced_ml <- bam(reduced_formula, data = dat, method = "ML")
  # Check that BMI/WHR dependence is not attributed to slope/WHR.
  bmi_whr <- bam(
    update(full_formula, . ~ . + ti(bmi, WHR)), data = dat, method = "fREML"
  )
  stopifnot(
    identical(full$model$value, reduced$model$value),
    nobs(full) == nrow(dat), all(is.finite(coef(full)))
  )
  smooth_test <- summary(full)$s.table["ti(Slope,WHR)", ]
  sensitivity_test <- summary(bmi_whr)$s.table["ti(Slope,WHR)", ]
  ml_test <- anova(reduced_ml, full_ml, test = "F")
  tests[[i]] <- data.frame(
    sex = sex_value, original_n = nrow(original$bmdf[[i]]), n = nrow(dat),
    whr_missing = nrow(original$bmdf[[i]]) - nrow(dat),
    whr_smooth_f = summary(full)$s.table["s(WHR)", "F"],
    whr_smooth_p = summary(full)$s.table["s(WHR)", "p-value"],
    interaction_edf = unname(smooth_test["edf"]),
    interaction_f = unname(smooth_test["F"]),
    interaction_p = unname(smooth_test["p-value"]),
    ml_deviance_improvement = ml_test[2, "Deviance"],
    ml_p = ml_test[2, "Pr(>F)"],
    bmi_whr_sensitivity_p = unname(sensitivity_test["p-value"]),
    r_squared_full = summary(full)$r.sq,
    r_squared_reduced = summary(reduced)$r.sq
  )
  diagnostics[[i]] <- data.frame(
    sex = sex_value, converged = isTRUE(full$converged),
    rank = full$rank, coefficients = length(coef(full)),
    min_k_check_p = min(k.check(full, n.rep = 400)[, "p-value"]),
    bmi_whr_correlation = cor(dat$bmi, dat$WHR),
    slope_whr_correlation = cor(dat$Slope, dat$WHR),
    mgcv_version = as.character(packageVersion("mgcv"))
  )
  fits[[i]] <- list(
    sex = sex_value, data = dat, full = full, reduced = reduced,
    full_ml = full_ml, reduced_ml = reduced_ml, bmi_whr = bmi_whr
  )

  # WHR profiles: sex-specific WHR quartiles among people near each BMI.
  sex_profiles <- bind_rows(lapply(bmi_values, function(bmi_value) {
    band <- filter(dat, abs(bmi - bmi_value) <= bmi_half_width)
    data.frame(
      sex = sex_value, bmi = bmi_value, age = profile_age,
      whr_profile = profile_levels,
      WHR = as.numeric(quantile(band$WHR, c(0.25, 0.5, 0.75))),
      bmi_band_n = nrow(band)
    )
  }))
  support <- bind_rows(lapply(seq_len(nrow(sex_profiles)), function(j) {
    local <- filter(
      dat, abs(bmi - sex_profiles$bmi[j]) <= bmi_half_width,
      abs(WHR - sex_profiles$WHR[j]) <= whr_half_width,
      abs(age - profile_age) <= age_half_width
    )
    data.frame(
      local_n = nrow(local),
      local_slope_p025 = as.numeric(quantile(local$Slope, 0.025)),
      local_slope_p975 = as.numeric(quantile(local$Slope, 0.975))
    )
  }))
  # Covariate values match the existing model-based contrasts, not a person.
  sex_profiles <- bind_cols(sex_profiles, support) |>
    mutate(
      CurrSmoke = 0.08, Diabetes = 0.05, Insulin = 0.5, LipidLower = 0.16,
      HTMed = 0.25, TSI = -1.3, CollegeEd = 0.3, Menopause = 0.4,
      ancestry = "EUR"
    )
  profiles[[i]] <- sex_profiles

  pred_grid <- crossing(sex_profiles, Slope = seq(-0.5, 0.5, by = 0.05)) |>
    mutate(
      within_local_slope_range =
        Slope >= local_slope_p025 & Slope <= local_slope_p975
    )
  predictions[[i]] <- bind_cols(
    pred_grid, log_liver_predictions(full, pred_grid)
  )

  # Upper versus lower slope quartile within each WHR profile.
  high <- low <- sex_profiles
  high$Slope <- pair[2]
  low$Slope <- pair[1]
  slope_contrasts[[i]] <- contrast_rows(
    full, high, low,
    select(sex_profiles, sex, bmi, age, whr_profile, WHR, local_n) |>
      mutate(
        slope_p25 = pair[1], slope_p75 = pair[2],
        within_local_slope_range = pair[1] >= sex_profiles$local_slope_p025 &
          pair[2] <= sex_profiles$local_slope_p975
      )
  )

  # Higher versus lower WHR at the same BMI and the same slope.
  slope_points <- c(p25 = pair[1], p50 = slope_median, p75 = pair[2])
  for (bmi_value in bmi_values) {
    hi <- filter(sex_profiles, bmi == bmi_value, whr_profile == "Higher (p75)")
    lo <- filter(sex_profiles, bmi == bmi_value, whr_profile == "Lower (p25)")
    support_low <- max(lo$local_slope_p025, hi$local_slope_p025)
    support_high <- min(lo$local_slope_p975, hi$local_slope_p975)
    for (nm in names(slope_points)) {
      hi$Slope <- lo$Slope <- slope_points[[nm]]
      profile_contrasts[[length(profile_contrasts) + 1]] <- contrast_rows(
        full, hi, lo,
        data.frame(
          sex = sex_value, bmi = bmi_value, age = profile_age,
          slope_percentile = nm, Slope = slope_points[[nm]],
          whr_p25 = lo$WHR, whr_p75 = hi$WHR,
          within_local_slope_range = slope_points[[nm]] >= support_low &
            slope_points[[nm]] <= support_high
        )
      )
    }
    # Ratio of ratios: slope-quartile contrast at higher versus lower WHR.
    differences[[length(differences) + 1]] <- bind_cols(
      data.frame(
        sex = sex_value, bmi = bmi_value, age = profile_age,
        whr_p25 = lo$WHR, whr_p75 = hi$WHR,
        slope_p25 = pair[1], slope_p75 = pair[2],
        within_local_slope_range = pair[1] >= support_low &
          pair[2] <= support_high
      ),
      log_liver_profile_difference(full, hi, lo, pair[1], pair[2])
    )
  }
}

tests <- bind_rows(tests) |>
  mutate(
    interaction_fdr = p.adjust(interaction_p, method = "BH"),
    ml_fdr = p.adjust(ml_p, method = "BH"),
    bmi_whr_sensitivity_fdr = p.adjust(bmi_whr_sensitivity_p, method = "BH")
  )
add_p <- function(x) mutate(x, p_value = 2 * pnorm(-abs(estimate / std_error)))
outputs <- list(
  tests = tests,
  diagnostics = bind_rows(diagnostics),
  profiles = bind_rows(profiles),
  predictions = bind_rows(predictions),
  slope_contrasts = add_p(bind_rows(slope_contrasts)),
  profile_contrasts = add_p(bind_rows(profile_contrasts)),
  profile_differences = add_p(bind_rows(differences))
)
stopifnot(
  nrow(outputs$slope_contrasts) == 18, nrow(outputs$profile_contrasts) == 18,
  nrow(outputs$profile_differences) == 6
)
for (nm in names(outputs)) {
  write_csv(
    outputs[[nm]], file.path(results_dir, paste0("liver_fat_whr_", nm, ".csv"))
  )
}
saveRDS(fits, model_path)
print(tests)
print(select(outputs$profile_contrasts, sex, bmi, slope_percentile, ratio,
             ratio_low, ratio_high))
print(select(outputs$profile_differences, sex, bmi, ratio_of_ratios,
             ratio_low, ratio_high))
# nolint end: object_usage_linter.

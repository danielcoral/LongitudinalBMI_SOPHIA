# Is the recent-slope association with liver fat independent of WHR?
# Reuses the saved script-35 liver-fat models that include a WHR smooth but
# no slope-by-WHR interaction, then fits the same model without WHR in the
# same participants. Compares the upper- versus lower-quartile slope contrast
# with and without WHR, and tests all slope terms once WHR is in the model.
# Inputs: liver_fat_whr_models.rds (processed_data), plus the script-35
#   grid and slope quartiles in liver_fat_whr_profiles.csv and
#   liver_fat_whr_slope_contrasts.csv (results).
# Output: liver_fat_slope_whr_adjustment.csv in data/UKB/results/.
# dplyr resolves column names through non-standard evaluation.
# nolint start: object_usage_linter.
library(mgcv)
library(readr)
library(dplyr, warn.conflicts = FALSE)
source("src/liver_fat_interaction_helpers.R")

results_dir <- "data/UKB/results"
fits <- readRDS("data/UKB/processed_data/liver_fat_whr_models.rds")
# Script-35 profiles: BMI 25/30/35 at age 60 with fixed covariate values.
# Without a slope-by-WHR term the slope contrast does not depend on WHR, so
# the median-WHR profile is used only to complete the prediction rows.
profiles <- read_csv(
  file.path(results_dir, "liver_fat_whr_profiles.csv"), show_col_types = FALSE
) |>
  filter(whr_profile == "Median (p50)")
quartiles <- read_csv(
  file.path(results_dir, "liver_fat_whr_slope_contrasts.csv"),
  show_col_types = FALSE
) |>
  distinct(sex, slope_p25, slope_p75)

rows <- list()
for (fit in fits) {
  with_whr <- fit$reduced
  labels <- attr(terms(formula(with_whr)), "term.labels")
  no_whr_formula <- reformulate(
    grep("WHR", labels, value = TRUE, invert = TRUE), response = "value"
  )
  no_slope_formula <- reformulate(
    grep("Slope", labels, value = TRUE, invert = TRUE), response = "value"
  )
  no_whr <- bam(no_whr_formula, data = fit$data, method = "fREML")
  no_slope_ml <- bam(no_slope_formula, data = fit$data, method = "ML")
  stopifnot(
    identical(no_whr$model$value, with_whr$model$value),
    !any(grepl("WHR", attr(terms(no_whr_formula), "term.labels")))
  )
  # Joint test of s(Slope), ti(bmi, Slope) and ti(Slope, age) given WHR.
  slope_test <- anova(no_slope_ml, fit$reduced_ml, test = "F")

  grid <- filter(profiles, sex == fit$sex)
  pair <- filter(quartiles, sex == fit$sex)
  high <- low <- grid
  high$Slope <- pair$slope_p75
  low$Slope <- pair$slope_p25
  without <- log_liver_contrast(no_whr, high, low)
  adjusted <- log_liver_contrast(with_whr, high, low)
  rows[[fit$sex]] <- data.frame(
    sex = fit$sex, n = nobs(with_whr), bmi = grid$bmi, age = grid$age,
    slope_p25 = pair$slope_p25, slope_p75 = pair$slope_p75,
    ratio_without_whr = without$ratio,
    ratio_without_whr_low = without$ratio_low,
    ratio_without_whr_high = without$ratio_high,
    ratio_with_whr = adjusted$ratio,
    ratio_with_whr_low = adjusted$ratio_low,
    ratio_with_whr_high = adjusted$ratio_high,
    # Share of the log-scale slope contrast removed by adjusting for WHR.
    attenuation_pct = 100 * (1 - adjusted$estimate / without$estimate),
    slope_deviance_given_whr = slope_test[2, "Deviance"],
    slope_p_given_whr = slope_test[2, "Pr(>F)"],
    r_squared_without_whr = summary(no_whr)$r.sq,
    r_squared_with_whr = summary(with_whr)$r.sq
  )
}

output <- bind_rows(rows) |>
  group_by(sex) |>
  mutate(slope_p_given_whr = first(slope_p_given_whr)) |>
  ungroup()
fdr <- distinct(output, sex, slope_p_given_whr) |>
  mutate(slope_fdr_given_whr = p.adjust(slope_p_given_whr, method = "BH"))
output <- left_join(output, fdr, by = c("sex", "slope_p_given_whr"))
stopifnot(nrow(output) == 6, all(is.finite(output$ratio_with_whr)))
write_csv(
  output, file.path(results_dir, "liver_fat_slope_whr_adjustment.csv")
)
print(select(
  output, sex, bmi, ratio_without_whr, ratio_with_whr,
  ratio_with_whr_low, ratio_with_whr_high, attenuation_pct,
  slope_p_given_whr, slope_fdr_given_whr
))
# nolint end: object_usage_linter.

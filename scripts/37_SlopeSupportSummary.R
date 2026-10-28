# Support of the within-cohort contrasts: age, BMI and slope distributions.
# Inputs: data/{ALSPAC,RS,UKB}/resultdf_*.rds (quantile summaries);
#   data/UKB/processed_data/lmmtabrec.rds and lmmtab2i.rds (individual UK
#   Biobank slopes; no individual-level ALSPAC or Rotterdam files are kept);
#   data/UKB/results/table_supp_08_cohorts.csv and
#   data/UKB/results/table_supp_09_slope_quantiles.csv (cross-checks).
# Output: data/cross_cohort/slope_support_summary.csv (one row per cohort,
#   sex and visit).
library(dplyr, warn.conflicts = FALSE)
library(purrr)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
results <- load_cohort_results(project_dir)

# bmsum rows 1-5 are the 2.5th, 25th, 50th, 75th and 97.5th percentiles.
support <- pmap_dfr(
  list(
    cohort = results$cohort, reftime = results$reftime, sex = results$sex,
    n = results$TotalN, bmsum = results$bmsum,
    long_summary = results$LongDatSum
  ),
  function(cohort, reftime, sex, n, bmsum, long_summary) {
    tibble(
      cohort, reftime, sex, n,
      followup_median = abs(long_summary$T_Md),
      age_p025 = reported_age(cohort, bmsum$age[1]),
      age_median = reported_age(cohort, bmsum$age[3]),
      age_p975 = reported_age(cohort, bmsum$age[5]),
      bmi_p025 = bmsum$bmi[1], bmi_p25 = bmsum$bmi[2],
      bmi_median = bmsum$bmi[3], bmi_p75 = bmsum$bmi[4],
      bmi_p975 = bmsum$bmi[5],
      slope_p025 = bmsum$Slope[1], slope_p25 = bmsum$Slope[2],
      slope_median = bmsum$Slope[3], slope_p75 = bmsum$Slope[4],
      slope_p975 = bmsum$Slope[5]
    )
  }
) |>
  mutate(
    slope_iqr = slope_p75 - slope_p25,
    zero_in_central95 = slope_p025 <= 0 & slope_p975 >= 0,
    minus05_in_central95 = slope_p025 <= -0.5 & slope_p975 >= -0.5,
    plus05_in_central95 = slope_p025 <= 0.5 & slope_p975 >= 0.5
  )

# Shares of individual UK Biobank slopes beyond +-0.5 kg/m2/year. The
# denominator n_slope is everyone with a slope estimate, so it can exceed the
# model-sample n slightly.
processed_dir <- file.path(project_dir, "data", "UKB", "processed_data")
shares <- imap_dfr(
  c(rec = "lmmtabrec.rds", `2i` = "lmmtab2i.rds"),
  function(file, visit) {
    models <- readRDS(file.path(processed_dir, file))
    map2_dfr(models$sex, models$indiv_est, function(sex, estimates) {
      tibble(
        cohort = "UKB", reftime = visit, sex, n_slope = nrow(estimates),
        share_below_minus05 = mean(estimates$Slope < -0.5),
        share_above_plus05 = mean(estimates$Slope > 0.5)
      )
    })
  }
)
support <- left_join(support, shares, by = c("cohort", "reftime", "sex"))

# Cross-check with the cohort tables (same source, so differences are errors).
results_dir <- file.path(project_dir, "data", "UKB", "results")
table_8 <- readr::read_csv(
  file.path(results_dir, "table_supp_08_cohorts.csv"), show_col_types = FALSE
)
table_9 <- readr::read_csv(
  file.path(results_dir, "table_supp_09_slope_quantiles.csv"),
  show_col_types = FALSE
)
keys <- c("cohort", "sex", "reftime")
check_8 <- inner_join(support, table_8, by = keys, suffix = c("", "_t8"))
check_9 <- inner_join(support, table_9, by = keys, suffix = c("", "_t9"))
stopifnot(
  nrow(support) == 8, nrow(check_8) == 8, nrow(check_9) == 8,
  all(check_8$n == check_8$n_t8),
  max(abs(check_8$bmi_median - check_8$bmi_median_t8)) < 1e-10,
  max(abs(check_9$slope_p25 - check_9$slope_p25_t9)) < 1e-10,
  max(abs(check_9$slope_p75 - check_9$slope_p75_t9)) < 1e-10
)

output_dir <- file.path(project_dir, "data", "cross_cohort")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(support, file.path(output_dir, "slope_support_summary.csv"))

# Extract manuscript numbers from existing model/data outputs, without refits.
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
library(mgcv)
library(lme4)

project_dir <- getwd()
processed_dir <- file.path(project_dir, "data", "UKB", "processed_data")
results_dir <- file.path(project_dir, "data", "UKB", "results")
source("src/biomarker_ranges.R")
source("src/within_cohort_helpers.R")
write_csv(biomarker_ranges, file.path(results_dir, "table_supp_10_ranges.csv"))

trajectory_summary <- map_dfr(
  c(primary = "lmmtabrec.rds", strict = "lmmtabrec_strict.rds",
    second_visit = "lmmtab1i.rds", imaging = "lmmtab2i.rds"),
  function(filename) {
    models <- readRDS(file.path(processed_dir, filename))
    per_person <- map2_dfr(models$LongD, models$sex, function(data, sex) {
      data |>
        group_by(eid) |>
        summarise(
          sex = sex, n_measures = n(), duration = abs(min(time)),
          .groups = "drop"
        )
    })
    bind_rows(
      per_person,
      mutate(per_person, sex = "Overall")
    ) |>
      group_by(sex) |>
      summarise(
        n = n(), median_measures = median(n_measures),
        percent_three_measures = 100 * mean(n_measures >= 3),
        followup_p25 = quantile(duration, 0.25),
        followup_median = median(duration),
        followup_p75 = quantile(duration, 0.75), .groups = "drop"
      ) |>
      left_join(map2_dfr(models$LongMod, models$sex, function(model, sex) {
        slope <- lme4::fixef(model)["time"]
        se <- sqrt(diag(vcov(model)))["time"]
        tibble(sex, mean_slope = slope, slope_ci_low = slope - 1.96 * se,
               slope_ci_high = slope + 1.96 * se)
      }), by = "sex")
  }, .id = "analysis"
)
write_csv(
  trajectory_summary, file.path(results_dir, "manuscript_trajectories.csv")
)

cohort_summary <- map_dfr(c("ALSPAC", "RS", "UKB"), function(cohort) {
  result <- readRDS(file.path(
    project_dir, "data", cohort, paste0("resultdf_", cohort, ".rds")
  ))
  if (!"reftime" %in% names(result)) result$reftime <- NA_character_
  pmap_dfr(
    list(
      sex = result$sex, reftime = result$reftime, total_n = result$TotalN,
      bmsum = result$bmsum, long_summary = result$LongDatSum
    ),
    function(sex, reftime, total_n, bmsum, long_summary) {
      tibble(
        cohort, sex, reftime, n = total_n,
        trajectory_n = long_summary$TotalInd,
        age_p25 = reported_age(cohort, bmsum$age[2]),
        age_median = reported_age(cohort, bmsum$age[3]),
        age_p75 = reported_age(cohort, bmsum$age[4]), bmi_p25 = bmsum$bmi[2],
        bmi_median = bmsum$bmi[3], bmi_p75 = bmsum$bmi[4],
        followup_median = abs(long_summary$T_Md),
        slope_p025 = bmsum$Slope[1], slope_p25 = bmsum$Slope[2],
        slope_median = bmsum$Slope[3], slope_p75 = bmsum$Slope[4],
        slope_p975 = bmsum$Slope[5]
      )
    }
  )
}) |>
  group_by(cohort, reftime) |>
  mutate(sex_percent = 100 * n / sum(n)) |>
  ungroup()
write_csv(cohort_summary, file.path(results_dir, "table_supp_08_cohorts.csv"))
write_csv(
  select(cohort_summary, cohort, sex, reftime, starts_with("slope_")),
  file.path(results_dir, "table_supp_09_slope_quantiles.csv")
)

model_summary <- map_dfr(
  c(primary = "slopebmregdf.rds", imaging = "slopebmregdf_2i.rds"),
  function(filename) {
    models <- readRDS(file.path(processed_dir, filename))
    tibble(
      sex = models$sex, trait = models$trait,
      n = map_int(models$mod1, nobs),
      variance_explained = map_dbl(models$mod1, ~summary(.x)$r.sq)
    )
  }, .id = "analysis"
)
write_csv(model_summary, file.path(results_dir, "manuscript_model_fit.csv"))
imaging_averages <- model_summary |>
  filter(analysis == "imaging") |>
  mutate(group = ifelse(trait == "LiverF_perc", "Liver", "DXA")) |>
  group_by(group) |>
  summarise(mean_variance_percent = 100 * mean(variance_explained),
            .groups = "drop")
write_csv(
  imaging_averages, file.path(results_dir, "imaging_variance_averages.csv")
)

bootstrap_inventory <- map_dfr(
  c("bootstrap_rec", "bootstrap_2i"), function(dir) {
    files <- list.files(file.path(processed_dir, dir), full.names = TRUE)
    map_dfr(files, function(file) {
      data <- read_tsv(file, show_col_types = FALSE, progress = FALSE)
      tibble(directory = dir, file = basename(file), rows = nrow(data),
             iterations = n_distinct(data$B_ind),
             complete = setequal(unique(data$B_ind), 1:1000))
    })
  }
)
write_csv(
  bootstrap_inventory, file.path(results_dir, "bootstrap_inventory.csv")
)
stopifnot(all(bootstrap_inventory$complete))

# Compact summaries for qualitative contrast statements in the Results.
contrast_summary <- read_csv(
  file.path(results_dir, "table_primary_rapid_gain_vs_loss.csv"),
  show_col_types = FALSE
) |>
  group_by(sex, trait) |>
  summarise(
    minimum_difference = min(estimate), maximum_difference = max(estimate),
    positive_over_point_one_sd = sum(estimate > thresh),
    negative_over_point_one_sd = sum(estimate < -thresh),
    evaluated_strata = n(), .groups = "drop"
  )
write_csv(
  contrast_summary, file.path(results_dir, "primary_contrast_summary.csv")
)

# Delay among eligible participants with a later biomarker measurement.
recruitment <- read_csv(
  file.path(project_dir, "data", "UKB", "raw_data",
            "Weight_Recruitment_participant.csv"), show_col_types = FALSE
) |>
  transmute(eid, recruitment_date = as.Date(p53_i0))
second <- read_csv(
  file.path(project_dir, "data", "UKB", "raw_data",
            "Weight_Instance1_participant.csv"), show_col_types = FALSE
) |>
  transmute(eid, second_date = as.Date(p53_i1))
lagged_data <- read_tsv(file.path(processed_dir, "bm1i.tsv"),
                        show_col_types = FALSE)
trajectory_models <- readRDS(file.path(processed_dir, "lmmtabrec.rds"))
eligible <- bind_rows(trajectory_models$indiv_est) |>
  select(eid) |>
  inner_join(lagged_data, by = "eid") |>
  inner_join(recruitment, by = "eid") |>
  inner_join(second, by = "eid")
write_csv(
  summarise(eligible, n = n(),
            median_delay_years = median(
              as.numeric(second_date - recruitment_date) / 365.25
            )),
  file.path(results_dir, "lagged_biomarker_delay.csv")
)

# Survival model samples: slope-only Cox models use all participants with
# slopes; full models (and GAMs) are complete cases for covariates/biomarkers,
# so they differ from the outcome summary table (Supplementary Table 18).
survival <- readRDS(file.path(processed_dir, "survival_models_current.rds"))
survival_samples <- survival |>
  transmute(
    sex, outcome,
    n_slope_only = map_int(survcox_slopeonly, "n"),
    events_slope_only = map_int(survcox_slopeonly, "nevent"),
    n_full = map_int(survcox_full, "n"),
    events_full = map_int(survcox_full, "nevent")
  )
write_csv(
  survival_samples, file.path(results_dir, "survival_model_samples.csv")
)

# Audit document references, including missing outputs, without hiding gaps.
documents <- c("paper/manus.md", "paper/main_figs.md",
               "paper/supp_figs.md", "paper/supp_tables.md")
references <- map_dfr(documents, function(document) {
  text <- readLines(document)
  map2_dfr(text, seq_along(text), function(line, number) {
    paths <- regmatches(line, gregexpr(
      "(?:data|plots)/[A-Za-z0-9_./-]+\\.(?:csv|tsv|txt|png|rds)",
      line, perl = TRUE
    ))[[1]]
    tibble(document, line = number, output = paths,
           exists = file.exists(paths))
  })
})
write_csv(
  references, file.path(results_dir, "manuscript_output_references.csv")
)
print(trajectory_summary)
print(cohort_summary)
print(imaging_averages)
print(filter(references, !exists))

# Generate the supplementary workbook from analysis outputs, not manual cells.
library(readr)
library(writexl)
project_dir <- getwd()
source_paths <- c(
  "data/cross_cohort/slope_support_summary.csv",
  "data/cross_cohort/within_cohort_contrasts.csv",
  "data/UKB/results/table_primary_slope_added_value.csv",
  "data/UKB/results/table_primary_interaction_added_value.csv",
  "data/UKB/results/imaging_slope_added_value.csv",
  "data/UKB/results/imaging_interaction_added_value.csv",
  "data/UKB/results/liver_fat_quartile_contrasts.csv",
  "data/UKB/results/liver_fat_whr_tests.csv",
  "data/UKB/results/liver_fat_whr_slope_contrasts.csv",
  "data/UKB/results/liver_fat_whr_profile_contrasts.csv",
  "data/UKB/results/liver_fat_whr_profile_differences.csv",
  "data/cross_cohort/stable_contrasts_common_grid.csv",
  "data/cross_cohort/stable_consistency_ukb_rs.csv",
  "data/UKB/results/table_primary_gain_vs_stable.csv",
  "data/UKB/results/table_imaging_gain_vs_stable.csv",
  "data/UKB/results/table_log_liver_fat_gain_vs_stable.csv",
  "data/UKB/results/slope_asymmetry_summary.csv",
  "data/UKB/results/table_supp_08_cohorts.csv",
  "data/UKB/results/table_supp_09_slope_quantiles.csv",
  "data/UKB/results/table_supp_10_ranges.csv",
  "data/UKB/results/participant_characteristics.csv",
  "data/UKB/results/trajectory_followup.csv",
  "data/UKB/results/outcome_summary.csv"
)
files <- file.path(project_dir, source_paths)
tables <- lapply(files, function(file) {
  read_csv(file, show_col_types = FALSE)
})
# Supplementary Table 19 reports central percentiles, not Tukey hinges.
names(tables) <- paste0("Supplementary Table ", seq_along(tables))
write_xlsx(tables, "paper/supp_tables_generated.xlsx")
writeLines(
  c("# Generated supplementary tables", "",
    "Workbook: supp_tables_generated.xlsx", "",
    paste0("Table ", seq_along(files), ": ", source_paths)),
  "paper/supp_tables_generated.md"
)

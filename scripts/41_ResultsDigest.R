# Results digest for the life-stage restructure: key numbers and the status of
# each planned claim (C1 to C9), computed from the outputs of scripts 36-40.
# Inputs: within_cohort_contrasts.csv, slope_support_summary.csv and
#   stable_consistency_ukb_rs.csv in data/cross_cohort/;
#   liver_fat_quartile_contrasts.csv, liver_fat_whr_profile_contrasts.csv,
#   liver_fat_whr_tests.csv and liver_fat_whr_profile_differences.csv
#   (all in data/UKB/results/).
# Output: results_digest.txt in data/cross_cohort/ (the text is also printed).
library(dplyr, warn.conflicts = FALSE)

project_dir <- getwd()
read_table <- function(folder, name) {
  readr::read_csv(
    file.path(project_dir, "data", folder, name), show_col_types = FALSE
  )
}
within <- read_table("cross_cohort", "within_cohort_contrasts.csv")
support <- read_table("cross_cohort", "slope_support_summary.csv")
agreement <- read_table("cross_cohort", "stable_consistency_ukb_rs.csv")
liver <- read_table("UKB/results", "liver_fat_quartile_contrasts.csv") |>
  filter(grid == "cohort_percentile")
tests <- read_table("UKB/results", "liver_fat_whr_tests.csv")
whr_profiles <- read_table(
  "UKB/results", "liver_fat_whr_profile_differences.csv"
)
whr_levels <- read_table("UKB/results", "liver_fat_whr_profile_contrasts.csv")
biomarkers <- c("ALT", "AST", "GGT", "TG", "LDL", "HDL", "SBP", "DBP",
                "GLU", "CRP", "SCR")
dxa <- c("BF_perc", "AndF_perc", "GynF_perc", "ArmsF_perc", "LegsF_perc",
         "TrunkF_perc")
verdict <- function(ok, marginal = FALSE) {
  if (ok) "SUPPORTED" else if (marginal) "MARGINAL" else "NOT SUPPORTED"
}

# Slope support: what the lower and upper quartile mean in each cohort.
slopes <- support |>
  transmute(
    cohort, visit = coalesce(reftime, "-"), sex, n, p2.5 = slope_p025,
    p25 = slope_p25, p50 = slope_median, p75 = slope_p75,
    p97.5 = slope_p975, iqr = slope_iqr
  ) |>
  mutate(across(p2.5:iqr, \(x) round(x, 2)))

# C1: direction of the contrast in 18 cells per trait (3 cohorts x 2 sexes x
# 3 BMI percentiles).
direction <- within |>
  filter(trait %in% biomarkers) |>
  group_by(trait) |>
  summarise(
    cells = n(), positive = sum(estimate_sd > 0),
    negative = sum(estimate_sd < 0),
    ci_excludes_0 = sum(conf_low_sd > 0 | conf_high_sd < 0)
  ) |>
  arrange(match(trait, biomarkers))
same_sign <- direction$trait[
  direction$positive == direction$cells |
    direction$negative == direction$cells
]
# C2: body-composition (DXA) contrasts, 18 cells per cohort and sex.
dxa_summary <- within |>
  filter(trait %in% dxa) |>
  group_by(cohort, sex) |>
  summarise(
    cells = n(), negative = sum(estimate_sd < 0),
    ci_excludes_0 = sum(conf_low_sd > 0 | conf_high_sd < 0),
    min_sd = round(min(estimate_sd), 3), max_sd = round(max(estimate_sd), 3),
    .groups = "drop"
  )
alspac_men <- filter(dxa_summary, cohort == "ALSPAC", sex == "Male")
other_dxa <- filter(dxa_summary, cohort != "ALSPAC")
# C4: liver-fat ratio should rise with BMI percentile within each sex.
rising <- liver |>
  arrange(sex, bmi) |>
  group_by(sex) |>
  summarise(rising = all(diff(estimate) > 0))
# C5: higher versus lower imaging-visit WHR at the same BMI and slope.
n_larger <- sum(whr_levels$ratio_low > 1)
# C6: interaction between slope and WHR.
male_p <- tests$interaction_p[tests$sex == "Male"]
female_p <- tests$interaction_p[tests$sex == "Female"]
# C7 and C8: where zero and +-0.5 kg/m2/year sit in each slope distribution.
alspac <- filter(support, cohort == "ALSPAC")
others <- filter(support, cohort != "ALSPAC")
alspac_ok <- all(alspac$slope_p25 > 0.3) &&
  !alspac$zero_in_central95[alspac$sex == "Male"] &&
  !any(alspac$minus05_in_central95)
others_ok <- all(
  others$zero_in_central95 & others$minus05_in_central95 &
    others$plus05_in_central95
)
# C9: same-sign agreement between UK Biobank and Rotterdam.
wording <- agreement |>
  filter(subset == "all_pairs") |>
  mutate(
    verdict = case_when(
      share_same_sign >= 0.8 ~ "consistently",
      share_same_sign >= 0.7 ~ "mostly",
      TRUE ~ "NOT consistent"
    )
  )

digest <- capture.output({
  cat("RESULTS DIGEST (scripts 36-40); contrasts are upper minus lower",
      "quartile of slope, in SD units unless stated\n")
  cat("\n[S] Slope distribution (kg/m2/year)\n")
  print(as.data.frame(slopes), row.names = FALSE)
  cat("\n[C1] Non-imaging traits: cells (of 18) with each sign\n")
  print(as.data.frame(direction), row.names = FALSE)
  cat("\n[C2] Body composition (DXA), 6 traits x 3 BMI percentiles\n")
  print(as.data.frame(dxa_summary), row.names = FALSE)
  cat("\n[C4] Liver fat ratio, upper/lower slope quartile (UK Biobank)\n")
  liver |>
    transmute(
      sex, pct = bmi_percentile, bmi = round(bmi, 1),
      ratio = round(ratio, 3), low = round(ratio_low, 3),
      high = round(ratio_high, 3), sd = round(estimate_sd, 3)
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
  cat("\n[C5] Liver-fat ratio, higher versus lower imaging-visit WHR\n")
  whr_levels |>
    transmute(
      sex, bmi, slope_percentile, ratio, low = ratio_low, high = ratio_high
    ) |>
    mutate(across(ratio:high, \(x) round(x, 3))) |>
    as.data.frame() |>
    print(row.names = FALSE)
  cat("\n[C6] Ratio of ratios, higher versus lower WHR profile\n")
  whr_profiles |>
    transmute(
      sex, bmi, ratio = round(ratio_of_ratios, 3), low = round(ratio_low, 3),
      high = round(ratio_high, 3)
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
  cat("\nCLAIM STATUS (use the fallback wording unless SUPPORTED)\n")
  cat("C1", verdict(length(same_sign) >= 3), "- same sign in all 18 cells:",
      paste(same_sign, collapse = ", "), "\n")
  cat("C2", verdict(
    alspac_men$negative == alspac_men$cells && all(other_dxa$negative == 0)
  ), sprintf(
    "- ALSPAC men %d of %d negative (%.3f to %.3f SD, %d CI excludes 0);",
    alspac_men$negative, alspac_men$cells, alspac_men$min_sd,
    alspac_men$max_sd, alspac_men$ci_excludes_0
  ), "RS and UKB all positive\n")
  cat("C3", verdict(all(support$slope_iqr >= 0.3 & support$slope_iqr <= 0.5)),
      "- slope IQR", round(min(support$slope_iqr), 2), "to",
      round(max(support$slope_iqr), 2), "kg/m2/year in all strata\n")
  cat("C4", verdict(all(liver$ratio > 1) && all(rising$rising)),
      "- ratio", round(min(liver$ratio), 3), "to", round(max(liver$ratio), 3),
      "rising with BMI in both sexes\n")
  cat("C5", verdict(n_larger == nrow(whr_levels)),
      "- higher WHR corresponds to higher liver fat (CI > 1) in", n_larger,
      "of", nrow(whr_levels), "profiles\n")
  cat("C6", verdict(min(male_p, female_p) < 0.05, min(male_p, female_p) < 0.1),
      "- slope x WHR interaction p: women", round(female_p, 3),
      "men", round(male_p, 3), "\n")
  cat("C7", verdict(alspac_ok),
      "- ALSPAC: lower quartile >0.3, -0.5 outside the central 95%, zero",
      "outside in men\n")
  cat("C8", verdict(others_ok),
      "- UK Biobank and Rotterdam: 0 and +-0.5 inside the central 95%\n")
  cat("C9", paste(
    sprintf("%s %d/%d %s", wording$comparison, wording$n_same_sign,
            wording$n_pairs, wording$verdict),
    collapse = "; "
  ), "\n")
})
writeLines(
  digest, file.path(project_dir, "data", "cross_cohort", "results_digest.txt")
)
writeLines(digest)

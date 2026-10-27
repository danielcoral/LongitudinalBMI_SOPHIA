# Directional agreement of the gain/stable/loss contrasts between UK Biobank
# and the Rotterdam Study on the common grid (age 65; BMI 22.5, 27.5, 32.5).
# One pair = one trait x sex x BMI x comparison present in both cohorts (17
# traits). "Same sign" compares the signs of the two estimates; the
# estimates share a model structure and are correlated, so only counts are
# reported, no p-values.
# Input: data/cross_cohort/stable_contrasts_common_grid.csv (script 40).
# Output: data/cross_cohort/stable_consistency_ukb_rs.csv (counts).
library(dplyr, warn.conflicts = FALSE)

project_dir <- getwd()
common <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "stable_contrasts_common_grid.csv"),
  show_col_types = FALSE
)
keys <- c("sex", "trait", "bmi", "comparison")
pairs <- inner_join(
  filter(common, cohort == "UKB"), filter(common, cohort == "RS"),
  by = keys, suffix = c("_ukb", "_rs")
) |>
  mutate(
    same_sign = sign(estimate_ukb) == sign(estimate_rs),
    ukb_excludes_zero = conf_low_ukb > 0 | conf_high_ukb < 0
  )
stopifnot(nrow(pairs) == 17 * 2 * 3 * 4, !anyNA(pairs$same_sign))

# Each pair is counted once per BMI value and once in the "all" row; the
# second subset keeps pairs whose UK Biobank interval excludes zero.
subsets <- bind_rows(
  mutate(pairs, subset = "all_pairs"),
  mutate(filter(pairs, ukb_excludes_zero), subset = "ukb_ci_excludes_zero")
)
agreement <- bind_rows(
  mutate(subsets, bmi = as.character(bmi)), mutate(subsets, bmi = "all")
) |>
  group_by(subset, comparison, bmi) |>
  summarise(
    n_pairs = n(), n_same_sign = sum(same_sign),
    share_same_sign = mean(same_sign), .groups = "drop"
  )

output_dir <- file.path(project_dir, "data", "cross_cohort")
readr::write_csv(
  agreement, file.path(output_dir, "stable_consistency_ukb_rs.csv")
)

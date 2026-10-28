# Directional agreement of the gain/stable/loss contrasts between UK Biobank
# (age 60) and the Rotterdam Study (age 70), both at BMI 30 kg/m2.
# One pair = one trait x sex x comparison present in both cohorts (17
# traits). "Same sign" compares the signs of the two estimates; only counts
# are reported, no p-values.
# Input: data/cross_cohort/stable_contrasts_common_grid.csv (script 39).
# Output: data/cross_cohort/stable_consistency_ukb_rs.csv (counts).
library(dplyr, warn.conflicts = FALSE)

project_dir <- getwd()
common <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "stable_contrasts_common_grid.csv"),
  show_col_types = FALSE
)
pairs <- inner_join(
  filter(common, cohort == "UKB"), filter(common, cohort == "RS"),
  by = c("sex", "trait", "comparison"), suffix = c("_ukb", "_rs")
) |>
  mutate(
    same_sign = sign(estimate_ukb) == sign(estimate_rs),
    ukb_excludes_zero = conf_low_ukb > 0 | conf_high_ukb < 0
  )
stopifnot(nrow(pairs) == 17 * 2 * 4, !anyNA(pairs$same_sign))

# The second subset keeps pairs whose UK Biobank interval excludes zero.
agreement <- bind_rows(
  mutate(pairs, subset = "all_pairs"),
  mutate(filter(pairs, ukb_excludes_zero), subset = "ukb_ci_excludes_zero")
) |>
  group_by(subset, comparison) |>
  summarise(
    n_pairs = n(), n_same_sign = sum(same_sign),
    share_same_sign = mean(same_sign), .groups = "drop"
  )

output_dir <- file.path(project_dir, "data", "cross_cohort")
readr::write_csv(
  agreement, file.path(output_dir, "stable_consistency_ukb_rs.csv")
)

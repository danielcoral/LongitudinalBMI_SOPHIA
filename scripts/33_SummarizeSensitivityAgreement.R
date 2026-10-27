# Export paired sensitivity agreement; avoid judging similarity from images.
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(purrr)
results_dir <- file.path(getwd(), "data", "UKB", "results")
settings <- tibble(
  file = c("table_sensitivity_excluding_selected_diagnoses.csv",
           "table_strict_followup_biomarker_comparison.csv",
           "table_second_visit_biomarker_comparison.csv"),
  column = c("Subset", "sample", "anchor"),
  reference = c("all", "all", "recruitment"),
  alternative = c("nodx", "strict", "secondvisit")
)
summary <- pmap_dfr(settings, function(file, column, reference, alternative) {
  data <- read_csv(file.path(results_dir, file), show_col_types = FALSE) |>
    transmute(sex, trait, age, bmi, subset = .data[[column]],
              estimate, conf.low, conf.high)
  paired <- data |>
    filter(subset == reference) |>
    select(-subset) |>
    inner_join(
      data |> filter(subset == alternative) |> select(-subset),
      by = c("sex", "trait", "age", "bmi"), suffix = c("_ref", "_alt")
    ) |>
    mutate(
      same_direction = sign(estimate_ref) == sign(estimate_alt),
      intervals_overlap = conf.low_ref <= conf.high_alt &
        conf.low_alt <= conf.high_ref,
      difference = estimate_alt - estimate_ref
    )
  write_csv(paired, file.path(results_dir, paste0(
    "sensitivity_paired_", alternative, ".csv"
  )))
  paired |>
    group_by(sex, trait) |>
    summarise(
      strata = n(), same_direction = sum(same_direction),
      intervals_overlap = sum(intervals_overlap),
      mean_absolute_difference = mean(abs(difference)),
      .groups = "drop"
    ) |>
    mutate(analysis = alternative, .before = 1)
})
write_csv(summary, file.path(results_dir, "sensitivity_agreement_summary.csv"))
print(summary)

# Summarize separate WHR models and the saved liver-fat/WHR interaction.
# Inputs: whr_trajectory_models.rds and liver_fat_whr_tests.csv.
# Outputs: whr_trajectory_summary.csv, whr_trajectory_contrasts.csv and
#   liver_fat_whr_added_value.csv in data/UKB/results/. No model refits.
library(mgcv)
library(readr)
library(dplyr, warn.conflicts = FALSE)
source("src/within_cohort_helpers.R")

results_dir <- "data/UKB/results"
fits <- readRDS("data/UKB/processed_data/whr_trajectory_models.rds")
summaries <- contrasts <- list()
for (i in seq_along(fits)) {
  fit <- fits[[i]]
  model <- fit$full
  dat <- model$model
  full_deviance <- deviance(model)
  explained <- model$null.deviance - full_deviance
  slope_test <- anova(fit$no_slope, model, test = "F")
  interaction_test <- anova(fit$no_interaction, model, test = "F")
  slope_p <- slope_test[2, 6]
  interaction_p <- interaction_test[2, 6]
  if (deviance(fit$no_slope) <= full_deviance) slope_p <- 1
  if (deviance(fit$no_interaction) <= full_deviance) interaction_p <- 1
  summaries[[i]] <- data.frame(
    visit = fit$visit, sex = fit$sex, n = nobs(model),
    variance_explained_pct = 100 * summary(model)$dev.expl,
    slope_p = slope_p,
    slope_share_explained_pct =
      100 * (deviance(fit$no_slope) - full_deviance) / explained,
    slope_added_variance_pct =
      100 * (deviance(fit$no_slope) - full_deviance) / model$null.deviance,
    interaction_p = interaction_p,
    age = median(dat$age),
    bmi_p25 = as.numeric(quantile(dat$bmi, 0.25)),
    bmi_p50 = median(dat$bmi),
    bmi_p75 = as.numeric(quantile(dat$bmi, 0.75)),
    slope_p25 = as.numeric(quantile(dat$Slope, 0.25)),
    slope_p75 = as.numeric(quantile(dat$Slope, 0.75))
  )
  row <- summaries[[i]]
  grid <- model_grid(
    model, row$age, c(row$bmi_p25, row$bmi_p50, row$bmi_p75)
  )
  estimates <- slope_pair_contrast(
    model, grid, row$slope_p25, row$slope_p75
  )
  contrasts[[i]] <- bind_cols(
    data.frame(
      visit = fit$visit, sex = fit$sex, n = row$n, age = row$age,
      bmi_percentile = c("p25", "p50", "p75"), bmi = grid$bmi,
      slope_p25 = row$slope_p25, slope_p75 = row$slope_p75
    ),
    estimates
  ) |>
    mutate(
      outcome_sd = outcome_sd(model),
      estimate_sd = estimate / outcome_sd,
      conf_low_sd = conf_low / outcome_sd,
      conf_high_sd = conf_high / outcome_sd
    )
}
summary_rows <- bind_rows(summaries) |>
  mutate(
    slope_fdr = p.adjust(slope_p, method = "BH"),
    interaction_fdr = p.adjust(interaction_p, method = "BH")
  )
contrast_rows <- bind_rows(contrasts) |>
  mutate(p_fdr = p.adjust(p_value, method = "BH"))
interaction <- read_csv(
  file.path(results_dir, "liver_fat_whr_tests.csv"),
  show_col_types = FALSE
) |>
  mutate(
    interaction_added_r_squared_pp =
      100 * (r_squared_full - r_squared_reduced)
  )
write_csv(
  summary_rows, file.path(results_dir, "whr_trajectory_summary.csv")
)
write_csv(
  contrast_rows, file.path(results_dir, "whr_trajectory_contrasts.csv")
)
write_csv(
  interaction, file.path(results_dir, "liver_fat_whr_added_value.csv")
)
print(summary_rows)
print(select(
  contrast_rows, visit, sex, bmi_percentile, estimate, conf_low, conf_high
))
print(select(
  interaction, sex, n, interaction_added_r_squared_pp, interaction_fdr
))

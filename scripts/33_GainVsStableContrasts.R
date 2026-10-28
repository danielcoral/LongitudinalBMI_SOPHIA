# Additional contrasts from saved models only: no refits or old output edits.
library(mgcv)
library(readr)
library(dplyr, warn.conflicts = FALSE)
source("src/slope_contrast_helpers.R")

processed_dir <- "data/UKB/processed_data"
results_dir <- "data/UKB/results"
model_files <- c(
  primary = "slopebmregdf.rds", imaging = "slopebmregdf_2i.rds",
  log_liver_fat = "slopebmregdf_2ilf.rds"
)
legacy_files <- file.path(
  results_dir, paste0("table_", names(model_files), "_rapid_gain_vs_loss.csv")
)
legacy_hashes <- tools::md5sum(legacy_files)
all_results <- list()
for (family in names(model_files)) {
  models <- readRDS(file.path(processed_dir, model_files[[family]]))
  rows <- lapply(seq_len(nrow(models)), function(i) {
    dat <- models$bmdf[[i]]
    result <- paired_slope_contrasts(models$mod1[[i]])
    # Preserve the original API's numerical-delta intervals exactly. The
    # analytic joint Vp calculation supplies the second-difference test.
    levels <- list(
      gain_vs_stable = c(0, 0.5), stable_vs_loss = c(-0.5, 0),
      gain_vs_loss = c(-0.5, 0.5)
    )
    grid <- slope_contrast_grid()
    grid$Slope <- NA_real_
    for (label in names(levels)) {
      api <- as.data.frame(marginaleffects::comparisons(
        models$mod1[[i]], newdata = grid,
        variables = list(Slope = levels[[label]])
      ))
      index <- which(result$comparison == label)
      stopifnot(
        identical(result$age[index], api$age),
        identical(result$bmi[index], api$bmi),
        max(abs(result$estimate[index] - api$estimate)) < 1e-8,
        max(abs(result$std_error[index] - api$std.error)) < 1e-4
      )
      result$std_error[index] <- api$std.error
      result$conf_low[index] <- api$conf.low
      result$conf_high[index] <- api$conf.high
      result$p_value[index] <- api$p.value
    }
    result <- mutate(
      result, model_family = family, sex = models$sex[i],
      trait = models$trait[i], n = nrow(dat), thresh = 0.1 * sd(dat$value),
      outcome_scale = ifelse(trait == "logLiverF_perc", "natural_log", "raw")
    )
    result |>
      mutate(
        ratio = ifelse(outcome_scale == "natural_log", exp(estimate), NA),
        ratio_low = ifelse(outcome_scale == "natural_log", exp(conf_low), NA),
        ratio_high = ifelse(outcome_scale == "natural_log", exp(conf_high), NA)
      )
  })
  results <- bind_rows(rows)
  legacy <- read_csv(
    file.path(results_dir, paste0("table_", family, "_rapid_gain_vs_loss.csv")),
    show_col_types = FALSE
  ) |>
    select(sex, trait, age, bmi, legacy_estimate = estimate,
           legacy_low = conf.low, legacy_high = conf.high)
  reproduction <- filter(results, comparison == "gain_vs_loss") |>
    inner_join(legacy, by = c("sex", "trait", "age", "bmi"))
  stopifnot(
    nrow(reproduction) == nrow(legacy),
    max(abs(reproduction$estimate - reproduction$legacy_estimate)) < 1e-6,
    max(abs(reproduction$conf_low - reproduction$legacy_low)) < 1e-6,
    max(abs(reproduction$conf_high - reproduction$legacy_high)) < 1e-6
  )
  write_csv(
    filter(results, comparison != "asymmetry"),
    file.path(results_dir, paste0("table_", family, "_gain_vs_stable.csv"))
  )
  all_results[[family]] <- results
  cat(family, ": legacy estimates and intervals reproduced.\n")
  rm(models)
  gc()
}
# The raw-liver model is copied into both imaging model tables.
combined <- bind_rows(all_results)
raw_liver <- filter(combined, trait == "LiverF_perc")
a <- filter(raw_liver, model_family == "imaging")
b <- filter(raw_liver, model_family == "log_liver_fat")
stopifnot(
  identical(a$estimate, b$estimate), identical(a$std_error, b$std_error),
  identical(a$conf_low, b$conf_low), identical(a$conf_high, b$conf_high)
)
# Keep duplicate raw-liver rows in reproduction tables, not in pooled counts.
combined <- filter(
  combined, !(model_family == "log_liver_fat" & trait == "LiverF_perc")
)
write_csv(
  combined, file.path(results_dir, "slope_reference_contrasts.csv")
)
stopifnot(identical(legacy_hashes, tools::md5sum(legacy_files)))
write_csv(
  data.frame(file = names(legacy_hashes), md5 = unname(legacy_hashes)),
  file.path(results_dir, "stable_reference_legacy_hashes.csv")
)
cat("All original contrast files unchanged. No model was refitted.\n")

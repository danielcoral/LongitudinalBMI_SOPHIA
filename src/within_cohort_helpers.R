# Within-cohort contrasts: upper versus lower quartile of recent BMI slope at
# the same BMI, computed from saved cohort models (no refitting).
# Needs harmonize_trait() from matched_cohort_helpers.R (source it first).
# dplyr resolves column names through non-standard evaluation.
# nolint start: object_usage_linter.
cohort_labels <- c(
  ALSPAC = "ALSPAC", UKB = "UK Biobank", RS = "Rotterdam Study"
)
cohort_colors <- c(ALSPAC = "#E69F00", UKB = "#0072B2", RS = "#009E73")
trait_groups <- tibble::tibble(
  trait = c(
    "ALT", "AST", "GGT", "HDL", "LDL", "TG", "SBP", "DBP", "GLU",
    "CRP", "SCR", "BF_perc", "AndF_perc", "GynF_perc",
    "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
  ),
  group = rep(
    c("Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney",
      "Imaging"),
    c(3, 3, 2, 1, 1, 1, 6)
  )
)
# Main traits measured in all three cohorts.
common_traits <- trait_groups$trait

# ALSPAC summaries use the earlier age; report the next assessment 5 years
# later. This changes reporting only, not the saved-model predictions.
reported_age <- function(cohort, age) {
  age + ifelse(cohort == "ALSPAC", 5, 0)
}

load_cohort_results <- function(project_dir) {
  cohorts <- c("ALSPAC", "RS", "UKB")
  files <- file.path(
    project_dir, "data", cohorts, paste0("resultdf_", cohorts, ".rds")
  )
  if (!all(file.exists(files))) {
    stop("Missing cohort result files: ",
         paste(files[!file.exists(files)], collapse = ", "))
  }
  purrr::map2_dfr(files, cohorts, function(path, cohort_name) {
    result <- readRDS(path)
    if (!"reftime" %in% names(result)) result$reftime <- NA_character_
    dplyr::mutate(result, cohort = cohort_name)
  })
}

cohort_models <- function(result_table) {
  models <- result_table |>
    tidyr::unnest(modres) |>
    dplyr::filter(
      !(cohort == "ALSPAC" & grepl("F_perc$", trait) &
          !grepl("TF_perc$", trait))
    ) |>
    dplyr::mutate(trait = harmonize_trait(trait)) |>
    # Cohort-specific extras (ALSPAC FI; Rotterdam Ca2, CysC, PO4, SHBG, TC,
    # Urate) and raw LiverF_perc are outside the comparison set; liver fat
    # is analysed on the log scale in script 38.
    dplyr::filter(trait %in% trait_groups$trait)
  if (anyDuplicated(models[c("cohort", "reftime", "sex", "trait")])) {
    stop("Duplicate cohort/visit/sex/trait models.")
  }
  models
}

# Covariates stay at each model's reference values (numeric: x[2], factor:
# first level), as in matched_cohort_contrasts().
model_grid <- function(model, age, bmi) {
  reference <- as.data.frame(lapply(
    model$var.summary,
    function(x) if (is.numeric(x)) x[2] else x[1]
  ))
  grid <- reference[rep(1, length(bmi)), , drop = FALSE]
  grid$age <- age
  grid$bmi <- bmi
  grid
}

# Paired difference in predictions between two slopes on one grid; the joint
# covariance Vp gives the interval.
slope_pair_contrast <- function(model, grid, slope_low, slope_high) {
  stopifnot(model$family$link == "identity")
  low <- high <- grid
  low$Slope <- slope_low
  high$Slope <- slope_high
  contrast <- predict(model, newdata = high, type = "lpmatrix") -
    predict(model, newdata = low, type = "lpmatrix")
  estimate <- as.vector(contrast %*% coef(model))
  std_error <- sqrt(rowSums((contrast %*% model$Vp) * contrast))
  data.frame(
    estimate = estimate, std_error = std_error,
    conf_low = estimate - qnorm(0.975) * std_error,
    conf_high = estimate + qnorm(0.975) * std_error,
    p_value = 2 * pnorm(-abs(estimate / std_error))
  )
}

# The saved models are bam fits with raw data stripped. For bam, df.null is
# the number of observations and null.deviance the total sum of squares, so
# the sample SD of the modelled outcome is recoverable.
outcome_sd <- function(model) {
  sqrt(model$null.deviance / (model$df.null - 1))
}

model_n <- function(model) {
  model$df.null
}

within_cohort_contrasts <- function(result_table) {
  models <- cohort_models(result_table)
  result <- purrr::pmap_dfr(
    models[c("cohort", "reftime", "sex", "trait", "mod1", "bmsum")],
    function(cohort, reftime, sex, trait, mod1, bmsum) {
      bmi <- as.numeric(bmsum$bmi[2:4])
      slope <- as.numeric(bmsum$Slope[c(2, 4)])
      age <- bmsum$age[3]
      if (any(!is.finite(c(bmi, slope, age))) || slope[2] <= slope[1]) {
        stop("Invalid BMI/slope/age summary for ", cohort, "/", sex)
      }
      cbind(
        tibble::tibble(
          cohort, reftime, sex, trait,
          bmi_percentile = c("p25", "p50", "p75"), bmi = bmi,
          age = reported_age(cohort, age),
          slope_p25 = slope[1], slope_p75 = slope[2]
        ),
        slope_pair_contrast(
          mod1, model_grid(mod1, age, bmi), slope[1], slope[2]
        ),
        outcome_sd = outcome_sd(mod1), n = model_n(mod1)
      )
    }
  )
  # ALSPAC android fat was supplied as lean percentage (100 - fat %). With an
  # identity link the fat % contrast is the negated estimate with a mirrored
  # interval. Do this before standardising by the outcome SD.
  inverted <- result$cohort == "ALSPAC" & result$trait == "AndF_perc"
  lower <- result$conf_low[inverted]
  result$estimate[inverted] <- -result$estimate[inverted]
  result$conf_low[inverted] <- -result$conf_high[inverted]
  result$conf_high[inverted] <- -lower
  result |>
    dplyr::group_by(cohort, sex, bmi_percentile) |>
    dplyr::mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      estimate_sd = estimate / outcome_sd,
      conf_low_sd = conf_low / outcome_sd,
      conf_high_sd = conf_high / outcome_sd,
      in_all_cohorts = trait %in% common_traits
    ) |>
    dplyr::left_join(trait_groups, by = "trait") |>
    dplyr::select(
      cohort, reftime, sex, trait, group, bmi_percentile, bmi, age,
      slope_p25, slope_p75, estimate, std_error, conf_low, conf_high,
      p_value, p_fdr, outcome_sd, estimate_sd, conf_low_sd, conf_high_sd,
      n, in_all_cohorts
    )
}
# nolint end: object_usage_linter.

# Fit separate WHR outcome models at recruitment and imaging.
# Inputs: processed anthropometry, covariates and saved BMI trajectories;
#   imaging-visit waist/hip measurements. Main analysis fits are unchanged.
# Saves model sets in data/UKB/processed_data/whr_trajectory_models.rds.
library(mgcv)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)

processed_dir <- "data/UKB/processed_data"
model_path <- file.path(processed_dir, "whr_trajectory_models.rds")
visits <- c("rec", "2i")
sexes <- c("Female", "Male")
whr_range <- c(0.6, 1.4)
covars <- c(
  "CurrSmoke", "Diabetes", "Insulin", "LipidLower", "HTMed",
  "TSI", "CollegeEd", "Menopause", "ancestry"
)
base_terms <- c(
  covars, "s(age, bs = 'cr')", "s(bmi, bs = 'cr')", "ti(age, bmi)"
)
slope_terms <- c("s(Slope, bs = 'cr')", "ti(Slope, age)")
full_formula <- reformulate(
  c(base_terms, slope_terms, "ti(bmi, Slope)"), response = "value"
)
reduced_formulas <- list(
  no_slope = reformulate(base_terms, response = "value"),
  no_interaction = reformulate(
    c(base_terms, slope_terms), response = "value"
  )
)
ancestry <- read_tsv(
  file.path(processed_dir, "ancestrytab.tsv"), col_types = "nf"
)
fits <- list()

for (visit in visits) {
  if (visit == "rec") {
    measurements <- read_tsv(
      file.path(processed_dir, "bmrec.tsv"), show_col_types = FALSE
    ) |>
      select(eid, sex, age, bmi, WHR)
  } else {
    waist_hip <- read_tsv(
      "data/UKB/raw_data/WaistHip_Instance2.tsv", show_col_types = FALSE
    ) |>
      transmute(eid = f.eid, WHR = f.48.2.0 / f.49.2.0)
    measurements <- read_tsv(
      file.path(processed_dir, "img2i.tsv"), show_col_types = FALSE
    ) |>
      select(eid, sex, age, bmi) |>
      inner_join(waist_hip, by = "eid")
  }
  cv <- read_tsv(
    file.path(processed_dir, paste0("cv", visit, ".tsv")),
    show_col_types = FALSE
  )
  trajectories <- readRDS(
    file.path(processed_dir, paste0("lmmtab", visit, ".rds"))
  ) |>
    select(sex, indiv_est) |>
    unnest(indiv_est)
  dat <- measurements |>
    inner_join(cv, by = "eid") |>
    inner_join(ancestry, by = "eid") |>
    inner_join(trajectories, by = c("eid", "sex")) |>
    select(eid, sex, age, bmi, Slope, all_of(covars), value = WHR) |>
    filter(is.finite(value), between(value, whr_range[1], whr_range[2])) |>
    drop_na()
  stopifnot(!anyDuplicated(dat$eid))

  for (sex_value in sexes) {
    sex_data <- dat |> filter(sex == sex_value) |> droplevels()
    cat("WHR:", visit, sex_value, "N =", nrow(sex_data), "\n")
    full <- bam(full_formula, data = sex_data, method = "fREML")
    reduced <- lapply(reduced_formulas, function(formula) {
      # Fixed full-model penalties match the main added-value analysis.
      setup <- bam(formula, data = sex_data, fit = FALSE)
      penalties <- full$sp[names(setup$sp)]
      stopifnot(all(is.finite(penalties)))
      bam(
        formula, data = sex_data, method = "fREML", sp = penalties
      )
    })
    stopifnot(
      isTRUE(full$converged),
      all(vapply(reduced, function(m) {
        isTRUE(m$converged) && identical(m$model$value, full$model$value)
      }, logical(1)))
    )
    fits[[paste(visit, sex_value, sep = "_")]] <- list(
      visit = visit, sex = sex_value, full = full,
      no_slope = reduced$no_slope,
      no_interaction = reduced$no_interaction
    )
  }
}
saveRDS(fits, model_path)
cat("Saved four WHR model sets to", model_path, "\n")

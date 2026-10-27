# Shared inclusion criteria; unbounded CRP lower limit matches legacy code.
biomarker_ranges <- tibble::tibble(
  trait = c("ALT", "AST", "CRP", "DBP", "GGT", "HDL", "LDL", "GLU",
            "SBP", "SCR", "TG"),
  lower = c(5, 15, -Inf, 50, 10, 0.5, 1, 3, 80, 40, 0.5),
  upper = c(70, 60, 10, 120, 120, 3, 7, 15, 220, 150, 7),
  unit = c("U/L", "U/L", "mg/L", "mmHg", "U/L", "mmol/L", "mmol/L",
           "mmol/L", "mmHg", "µmol/L", "mmol/L")
)
filter_biomarker_ranges <- function(data) {
  for (i in seq_len(nrow(biomarker_ranges))) {
    trait <- biomarker_ranges$trait[i]
    value <- data[[trait]]
    data[[trait]] <- ifelse(
      value >= biomarker_ranges$lower[i] &
        value <= biomarker_ranges$upper[i], value, NA_real_
    )
  }
  data
}

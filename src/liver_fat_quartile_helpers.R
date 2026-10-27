# Liver-fat contrasts between two slopes (the lower and upper quartile of
# recent BMI slope) on the log scale; ratio = exp(estimate) = upper / lower.
# log_liver_contrast() comes from liver_fat_interaction_helpers.R (source it
# first).
# nolint start: object_usage_linter.
quartile_contrast_rows <- function(model, grid, slope_pair) {
  grid$Slope <- NULL
  high <- low <- grid
  low$Slope <- slope_pair[1]
  high$Slope <- slope_pair[2]
  dplyr::bind_cols(
    dplyr::mutate(grid, slope_p25 = slope_pair[1], slope_p75 = slope_pair[2]),
    log_liver_contrast(model, high, low)
  )
}
# nolint end: object_usage_linter.

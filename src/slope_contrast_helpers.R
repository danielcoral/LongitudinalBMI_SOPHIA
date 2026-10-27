# Paired differences and curvature use one model's joint covariance.
slope_contrast_grid <- function() {
  expand.grid(
    CurrSmoke = 0.08, Diabetes = 0.05, LipidLower = 0.16,
    HTMed = 0.25, Insulin = 0.5, TSI = -1.3, CollegeEd = 0.3,
    Menopause = 0.4, age = c(45, 55, 65),
    bmi = c(22.5, 27.5, 32.5, 37.5), ancestry = "EUR",
    stringsAsFactors = FALSE
  )
}

paired_slope_contrasts <- function(model, grid = slope_contrast_grid()) {
  stopifnot(model$family$link == "identity")
  loss <- stable <- gain <- grid
  loss$Slope <- -0.5
  stable$Slope <- 0
  gain$Slope <- 0.5
  x_loss <- predict(model, newdata = loss, type = "lpmatrix")
  x_stable <- predict(model, newdata = stable, type = "lpmatrix")
  x_gain <- predict(model, newdata = gain, type = "lpmatrix")
  matrices <- list(
    gain_vs_stable = x_gain - x_stable,
    stable_vs_loss = x_stable - x_loss,
    gain_vs_loss = x_gain - x_loss,
    asymmetry = x_gain + x_loss - 2 * x_stable
  )
  estimates <- lapply(names(matrices), function(comparison) {
    contrast <- matrices[[comparison]]
    estimate <- as.vector(contrast %*% coef(model))
    std_error <- sqrt(pmax(0, rowSums((contrast %*% model$Vp) * contrast)))
    data.frame(
      grid, comparison = comparison, estimate = estimate,
      std_error = std_error,
      conf_low = estimate - qnorm(0.975) * std_error,
      conf_high = estimate + qnorm(0.975) * std_error,
      p_value = 2 * pnorm(-abs(estimate / std_error)),
      row.names = NULL
    )
  })
  do.call(rbind, estimates)
}

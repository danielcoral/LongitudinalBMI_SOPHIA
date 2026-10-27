# Shared prediction and paired-contrast calculations on the log-liver scale.
log_liver_predictions <- function(model, newdata) {
  design <- predict(model, newdata = newdata, type = "lpmatrix")
  estimate <- as.vector(design %*% coef(model))
  std_error <- sqrt(rowSums((design %*% model$Vp) * design))
  data.frame(
    log_estimate = estimate,
    log_std_error = std_error,
    liver_fat = exp(estimate),
    conf_low = exp(estimate - qnorm(0.975) * std_error),
    conf_high = exp(estimate + qnorm(0.975) * std_error)
  )
}

log_liver_contrast <- function(model, high, low) {
  contrast <- predict(model, newdata = high, type = "lpmatrix") -
    predict(model, newdata = low, type = "lpmatrix")
  estimate <- as.vector(contrast %*% coef(model))
  std_error <- sqrt(rowSums((contrast %*% model$Vp) * contrast))
  data.frame(
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - qnorm(0.975) * std_error,
    conf_high = estimate + qnorm(0.975) * std_error,
    ratio = exp(estimate),
    ratio_low = exp(estimate - qnorm(0.975) * std_error),
    ratio_high = exp(estimate + qnorm(0.975) * std_error)
  )
}

# Compare trajectory contrasts between body-fat profiles using joint Vp.
log_liver_profile_difference <- function(model, high_bf, low_bf,
                                         slope_low, slope_high) {
  high_gain <- high_loss <- high_bf
  low_gain <- low_loss <- low_bf
  high_gain$Slope <- low_gain$Slope <- slope_high
  high_loss$Slope <- low_loss$Slope <- slope_low
  contrast <- predict(model, high_gain, type = "lpmatrix") -
    predict(model, high_loss, type = "lpmatrix") -
    predict(model, low_gain, type = "lpmatrix") +
    predict(model, low_loss, type = "lpmatrix")
  estimate <- as.vector(contrast %*% coef(model))
  std_error <- sqrt(rowSums((contrast %*% model$Vp) * contrast))
  data.frame(
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - qnorm(0.975) * std_error,
    conf_high = estimate + qnorm(0.975) * std_error,
    ratio_of_ratios = exp(estimate),
    ratio_low = exp(estimate - qnorm(0.975) * std_error),
    ratio_high = exp(estimate + qnorm(0.975) * std_error)
  )
}

# Plot BMI/slope density regions and a regression line from saved trajectories.
# Inputs: sex-specific trajectory table and visit-specific BMI axis label.
# Output: ggplot object for the supplementary trajectory distributions.
# dplyr/ggplot resolve data columns through non-standard evaluation.
# nolint start: object_usage_linter.
trajectory_density_plot <- function(trajectories, bmi_label) {
  trajectories |>
    transmute(
      sex,
      plot_data = purrr::map(LongD, select, eid, bmi0),
      plot_data = purrr::map(plot_data, distinct),
      plot_data = purrr::map2(plot_data, indiv_est, inner_join, by = "eid")
    ) |>
    tidyr::unnest(plot_data) |>
    ggplot(aes(x = bmi0, y = Slope)) +
    ggdensity::geom_hdr() +
    geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
    facet_wrap(~sex) +
    theme_bw() +
    labs(
      x = bmi_label, y = "Individual slope (kg/m² per year)",
      alpha = "Probability\ndensity"
    ) +
    theme(
      legend.title = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      axis.text = element_text(size = 6)
    ) +
    coord_cartesian(xlim = c(20, 45), ylim = c(-2.5, 2.5))
}
# nolint end: object_usage_linter.

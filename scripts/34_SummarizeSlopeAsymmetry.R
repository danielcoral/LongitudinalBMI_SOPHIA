# Summarize signed curvature, retaining both increments and their uncertainty.
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(ggplot2)

results_dir <- "data/UKB/results"
contrasts <- read_csv(
  file.path(results_dir, "slope_reference_contrasts.csv"),
  show_col_types = FALSE
)
components <- contrasts |>
  filter(comparison %in% c("gain_vs_stable", "stable_vs_loss")) |>
  select(
    model_family, sex, trait, age, bmi, comparison,
    estimate, std_error, conf_low, conf_high
  ) |>
  pivot_wider(
    names_from = comparison,
    values_from = c(estimate, std_error, conf_low, conf_high)
  )
asymmetry <- contrasts |>
  filter(comparison == "asymmetry") |>
  left_join(components, by = c("model_family", "sex", "trait", "age", "bmi")) |>
  group_by(model_family) |>
  mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(
    gain_increment_larger = estimate_gain_vs_stable > estimate_stable_vs_loss,
    gain_ci_excludes_zero = conf_low_gain_vs_stable > 0 |
      conf_high_gain_vs_stable < 0,
    loss_ci_excludes_zero = conf_low_stable_vs_loss > 0 |
      conf_high_stable_vs_loss < 0,
    asymmetry_ci_excludes_zero = conf_low > 0 | conf_high < 0
  )
write_csv(asymmetry, file.path(results_dir, "slope_asymmetry_summary.csv"))

focus <- filter(
  asymmetry, trait %in% c("ALT", "AST", "GGT", "TG", "SBP", "logLiverF_perc"),
  bmi %in% c(32.5, 37.5)
)
write_csv(focus, file.path(results_dir, "slope_asymmetry_focus.csv"))
counts <- asymmetry |>
  group_by(model_family, sex, trait) |>
  summarise(
    n_strata = n(), gain_increment_larger = sum(gain_increment_larger),
    ci_excludes_zero = sum(asymmetry_ci_excludes_zero),
    fdr_below_005 = sum(p_fdr < 0.05), .groups = "drop"
  )
write_csv(counts, file.path(results_dir, "slope_asymmetry_counts.csv"))
writeLines(
  c(
    "Exploratory asymmetry = gain-minus-stable minus stable-minus-loss.",
    "This is signed curvature at slopes -0.5/0/+0.5, not causal reversal.",
    "Intervals use joint coefficient covariance, conditional on saved fits.",
    "BH correction: all sex/trait/age/BMI asymmetry tests within model family.",
    "Raw liver fat is counted once; log liver fat is a separate scale.",
    capture.output(print(as.data.frame(counts), row.names = FALSE)),
    "High-BMI focus (exact estimates and intervals also in CSV):",
    capture.output(print(as.data.frame(focus), row.names = FALSE))
  ),
  file.path(results_dir, "results_slope_asymmetry.txt")
)

plot_data <- contrasts |>
  filter(comparison %in% c("gain_vs_stable", "stable_vs_loss")) |>
  mutate(
    comparison = factor(
      comparison, c("gain_vs_stable", "stable_vs_loss"),
      c("Gain vs stable", "Stable vs loss")
    ),
    age = factor(age),
    estimate_sd = estimate / (10 * thresh),
    low_sd = conf_low / (10 * thresh), high_sd = conf_high / (10 * thresh)
  )
colors <- c("#0072B2", "#D55E00")
for (family in c("primary", "imaging", "log_liver_fat")) {
  values <- filter(plot_data, model_family == family)
  if (family == "log_liver_fat") {
    values <- filter(values, trait == "logLiverF_perc")
  }
  figure <- ggplot(
    values,
    aes(bmi, estimate_sd, color = comparison, shape = age,
        group = interaction(comparison, age))
  ) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.1, ymax = 0.1,
             fill = "grey93") +
    geom_hline(yintercept = 0, color = "grey60") +
    geom_errorbar(
      aes(ymin = low_sd, ymax = high_sd),
      position = position_dodge(width = 1.4), width = 0.3, linewidth = 0.3
    ) +
    geom_point(position = position_dodge(width = 1.4), size = 1.4) +
    facet_wrap(vars(trait, sex), ncol = 4) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = c(22.5, 27.5, 32.5, 37.5)) +
    labs(
      x = expression("Current BMI (kg/" * m^2 * ")"),
      y = "Model-based difference (sex-specific outcome SD)",
      color = NULL, shape = "Age",
      caption = paste(
        "Slopes: gain +0.5, stable 0, loss -0.5 kg/m²/year.\n",
        "Pointwise 95% CIs; shading +/-0.1 SD."
      )
    ) +
    theme_bw(base_size = 9.5) +
    theme(
      legend.position = "bottom", panel.grid.minor = element_blank(),
      axis.text = element_text(size = 9.5),
      strip.text = element_text(size = 9.5),
      legend.text = element_text(size = 9.5),
      legend.title = element_text(size = 9.5),
      plot.caption = element_text(size = 9.5, hjust = 0)
    )
  filename <- switch(
    family, primary = "gain_vs_stable.png",
    imaging = "gain_vs_stable_imaging.png",
    log_liver_fat = "gain_vs_stable_log_liver_fat.png"
  )
  ggsave(
    file.path("plots", filename), figure, width = 7,
    height = switch(family, primary = 8.4, imaging = 6.4, log_liver_fat = 3.5),
    units = "in", dpi = 300, bg = "white"
  )
}
print(counts)

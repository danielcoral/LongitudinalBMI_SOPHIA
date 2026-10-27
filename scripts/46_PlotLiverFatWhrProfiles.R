# Figure 4: liver fat by waist-to-hip ratio (WHR) in UK Biobank, at the same
# BMI and recent BMI slope (age 60; WHR and liver fat at the imaging visit).
# A: predicted liver fat across recent BMI slope for lower, median and higher
#    WHR profiles; dashed lines mark the lower and upper slope quartile.
# B: liver-fat ratio, higher (p75) versus lower (p25) WHR, at the same BMI and
#    the same slope (lower quartile, median, upper quartile).
# C: ratio of ratios: slope-quartile contrast at higher versus lower WHR.
# Inputs: liver_fat_whr_predictions.csv, liver_fat_whr_profile_contrasts.csv,
#   liver_fat_whr_profile_differences.csv and liver_fat_whr_slope_contrasts.csv
#   in data/UKB/results/ (script 36).
# Output figure: plots/liver_fat_whr_profiles.png
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)

source("src/figure_helpers.R")
read_result <- function(name) {
  readr::read_csv(file.path("data", "UKB", "results", name),
                  show_col_types = FALSE)
}
profile_levels <- c("Lower (p25)", "Median (p50)", "Higher (p75)")
profile_colors <- setNames(c("#0072B2", "#666666", "#D55E00"), profile_levels)

predictions <- read_result("liver_fat_whr_predictions.csv") |>
  mutate(
    whr_profile = factor(whr_profile, profile_levels),
    bmi_label = paste("BMI", bmi)
  )
quartiles <- read_result("liver_fat_whr_slope_contrasts.csv") |>
  distinct(sex, slope_p25, slope_p75) |>
  tidyr::pivot_longer(c(slope_p25, slope_p75), values_to = "slope")
profile_contrasts <- read_result("liver_fat_whr_profile_contrasts.csv") |>
  mutate(slope_percentile = factor(
    slope_percentile, c("p25", "p50", "p75"),
    c("Lower quartile", "Median", "Upper quartile")
  ))
differences <- read_result("liver_fat_whr_profile_differences.csv")

plot_a <- ggplot(predictions, aes(Slope, liver_fat, colour = whr_profile)) +
  geom_vline(
    data = quartiles, aes(xintercept = slope), linetype = "dashed",
    colour = "grey40", linewidth = 0.3
  ) +
  # Solid lines and ribbons: slopes inside the local range of the data.
  geom_ribbon(
    data = filter(predictions, within_local_slope_range),
    aes(ymin = conf_low, ymax = conf_high, fill = whr_profile),
    alpha = 0.1, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.6, linetype = "dotted") +
  geom_line(
    data = filter(predictions, within_local_slope_range), linewidth = 0.6
  ) +
  facet_grid(sex ~ bmi_label) +
  scale_colour_manual(values = profile_colors) +
  scale_fill_manual(values = profile_colors) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
  labs(
    x = expression("Recent BMI slope (kg/" * m^2 * "/year)"),
    y = "Geometric-mean\nliver fat (%)",
    colour = "WHR within the BMI range"
  ) +
  paper_theme

plot_b <- ggplot(
  profile_contrasts, aes(factor(bmi), ratio, colour = slope_percentile)
) +
  geom_hline(
    yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3
  ) +
  geom_pointrange(
    aes(ymin = ratio_low, ymax = ratio_high),
    position = position_dodge(width = 0.5), size = 0.25, linewidth = 0.4
  ) +
  facet_grid(. ~ sex) +
  scale_colour_manual(values = c("#56B4E9", "#666666", "#E69F00")) +
  labs(
    x = expression("BMI (kg/" * m^2 * ")"),
    y = "Liver-fat ratio\n(WHR p75 / p25)", colour = "Recent BMI slope"
  ) +
  paper_theme

plot_c <- ggplot(differences, aes(factor(bmi), ratio_of_ratios)) +
  geom_hline(
    yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3
  ) +
  geom_pointrange(
    aes(ymin = ratio_low, ymax = ratio_high), size = 0.25, linewidth = 0.4
  ) +
  facet_grid(. ~ sex) +
  labs(
    x = expression("BMI (kg/" * m^2 * ")"),
    y = "Ratio of slope-quartile\nratios (WHR p75 / p25)"
  ) +
  paper_theme

figure <- plot_a / (plot_b | plot_c) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure, "liver_fat_whr_profiles.png", 7, 7.5)

# Figure 3: liver fat by recent BMI slope.
# A: upper minus lower quartile of slope in SD units at each cohort's own BMI
#    percentiles, with ALT and body fat % as reference traits.
# B: liver-fat ratio (upper / lower quartile) on the fixed age x BMI grid,
#    supported cells only (age and BMI inside the central 95% of the sample).
# Inputs: data/cross_cohort/within_cohort_contrasts.csv (script 37);
#   data/UKB/results/liver_fat_quartile_contrasts.csv (script 39).
# Output figure: plots/liver_fat_quartile_contrasts.png
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "figure_helpers.R"))
liver <- readr::read_csv(
  file.path(project_dir, "data", "UKB", "results",
            "liver_fat_quartile_contrasts.csv"),
  show_col_types = FALSE
)
within <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "within_cohort_contrasts.csv"),
  show_col_types = FALSE
)
measures <- c("logLiverF_perc", "ALT", "BF_perc")
panel_a <- bind_rows(
  filter(within, cohort %in% c("UKB", "RS"), trait %in% measures),
  filter(liver, grid == "cohort_percentile")
) |>
  mutate(
    trait = factor(trait, rev(measures)),
    cohort = factor(cohort, c("UKB", "RS"))
  )
percentile_labels <- c(
  p25 = "BMI percentile 25", p50 = "BMI percentile 50",
  p75 = "BMI percentile 75"
)
plot_a <- ggplot(panel_a, aes(estimate_sd, trait, colour = cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(
    aes(xmin = conf_low_sd, xmax = conf_high_sd),
    position = position_dodge(width = 0.6), size = 0.2, linewidth = 0.35
  ) +
  facet_grid(sex ~ bmi_percentile, labeller = labeller(
    bmi_percentile = percentile_labels
  )) +
  scale_y_discrete(labels = trait_labels) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1, 0.2)) +
  coord_cartesian(xlim = c(-0.12, 0.33)) +
  scale_colour_manual(values = cohort_colors, labels = cohort_labels) +
  labs(
    x = "Upper minus lower quartile of slope (SD units)", y = NULL,
    colour = NULL
  ) +
  paper_theme

# Cells outside the central 95% of age or BMI are not drawn.
panel_b <- liver |>
  filter(grid == "ukb_fixed", age_in_central95, bmi_in_central95) |>
  mutate(age = factor(age))
plot_b <- ggplot(panel_b, aes(bmi, ratio, colour = age, group = age)) +
  geom_hline(
    yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3
  ) +
  geom_line(position = position_dodge(width = 1), linewidth = 0.3) +
  geom_pointrange(
    aes(ymin = ratio_low, ymax = ratio_high),
    position = position_dodge(width = 1), size = 0.2, linewidth = 0.35
  ) +
  facet_grid(. ~ sex, scales = "free_x", space = "free_x") +
  scale_colour_manual(
    values = c("45" = "#9ECAE1", "55" = "#4292C6", "65" = "#08519C")
  ) +
  scale_x_continuous(breaks = c(22.5, 27.5, 32.5, 37.5)) +
  labs(
    x = expression("BMI (kg/" * m^2 * ")"),
    y = "Liver-fat ratio\n(upper / lower slope quartile)",
    colour = "Age (years)"
  ) +
  paper_theme

figure <- plot_a / plot_b +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure, "liver_fat_quartile_contrasts.png", 7, 7)

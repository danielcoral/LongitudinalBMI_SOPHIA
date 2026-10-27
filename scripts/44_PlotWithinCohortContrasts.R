# Figure 2: lower versus upper quartile of recent BMI slope at the same BMI in
# ALSPAC, UK Biobank and the Rotterdam Study (the 17 traits measured in all
# three cohorts). Differences are in standard deviations (SD) of each trait.
# Input: data/cross_cohort/within_cohort_contrasts.csv (script 37).
# Output figure: plots/within_cohort_contrasts.png
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "figure_helpers.R"))
contrasts <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "within_cohort_contrasts.csv"),
  show_col_types = FALSE
) |>
  filter(in_all_cohorts) |>
  mutate(
    cohort = factor(cohort, names(cohort_colors)),
    trait = factor(trait, common_traits)
  )

# Block A: liver, lipids, blood pressure and glucose (9 traits). Block B is the
# same plot with the other 8 traits (inflammation, kidney, body composition).
block_a <- ggplot(
  filter(contrasts, trait %in% common_traits[1:9]),
  aes(bmi_percentile, estimate_sd, colour = cohort)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(
    aes(ymin = conf_low_sd, ymax = conf_high_sd),
    position = position_dodge(width = 0.7), size = 0.2, linewidth = 0.35
  ) +
  facet_grid(trait ~ sex, labeller = labeller(trait = trait_labels)) +
  scale_colour_manual(values = cohort_colors, labels = cohort_labels) +
  scale_x_discrete(labels = c("25", "50", "75")) +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25)) +
  coord_cartesian(ylim = c(-0.42, 0.5)) +
  labs(
    x = "BMI percentile within cohort and sex",
    y = "Upper minus lower quartile of slope (SD units)", colour = NULL
  ) +
  paper_theme +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5))
block_b <- block_a +
  filter(contrasts, trait %in% common_traits[10:17]) +
  labs(y = NULL)

figure <- (block_a | block_b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
save_figure(figure, "within_cohort_contrasts.png", 7, 7.2)

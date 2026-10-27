# Figure 1: the life-stage windows of the three cohorts and the distribution
# of the 5-year BMI slope in each. The thick bar is the interquartile range:
# its ends are the lower and upper quartile used for the contrasts.
# Input: data/cross_cohort/slope_support_summary.csv (script 38).
# Output figure: plots/design_slope_distributions.png
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "figure_helpers.R"))
strata <- c(
  "ALSPAC", "UK Biobank, recruitment", "UK Biobank, imaging visit",
  "Rotterdam Study"
)
support <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "slope_support_summary.csv"),
  show_col_types = FALSE
) |>
  mutate(
    stratum = case_when(
      cohort == "ALSPAC" ~ strata[1],
      cohort == "UKB" & reftime == "rec" ~ strata[2],
      cohort == "UKB" ~ strata[3],
      TRUE ~ strata[4]
    ),
    stratum = factor(stratum, rev(strata))
  )
cohort_scale <- scale_colour_manual(values = cohort_colors, guide = "none")

# A: ages covered (median and central 95%), the two sexes side by side.
dodge <- position_dodge(width = 0.55)
age_plot <- ggplot(support, aes(y = stratum, colour = cohort, shape = sex)) +
  geom_linerange(
    aes(xmin = age_p025, xmax = age_p975), position = dodge, linewidth = 0.5
  ) +
  geom_point(aes(x = age_median), position = dodge, size = 1.8) +
  cohort_scale +
  scale_shape_manual(values = c(Female = 16, Male = 17)) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  labs(
    x = "Age at biomarker assessment (years)", y = NULL, shape = NULL,
    title = "Age at assessment: median and central 95%"
  ) +
  paper_theme +
  theme(
    legend.position = "inside", legend.position.inside = c(0.8, 0.86),
    legend.background = element_blank(), legend.direction = "horizontal"
  )

# B: slope distribution; dotted lines mark -0.5 and +0.5 kg/m2/year.
slope_plot <- ggplot(support, aes(y = stratum, colour = cohort)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_vline(
    xintercept = c(-0.5, 0.5), colour = "grey50", linetype = "dotted",
    linewidth = 0.3
  ) +
  geom_linerange(aes(xmin = slope_p025, xmax = slope_p975), linewidth = 0.5) +
  geom_linerange(aes(xmin = slope_p25, xmax = slope_p75), linewidth = 2.6) +
  geom_point(aes(x = slope_median), shape = 21, fill = "white", size = 1.6) +
  geom_text(
    aes(x = slope_p25, label = sprintf("%.2f", slope_p25)),
    vjust = -1.1, hjust = 1, size = 2.5, colour = "grey20"
  ) +
  geom_text(
    aes(x = slope_p75, label = sprintf("%.2f", slope_p75)),
    vjust = -1.1, hjust = 0, size = 2.5, colour = "grey20"
  ) +
  facet_grid(sex ~ .) +
  cohort_scale +
  scale_x_continuous(breaks = seq(-1, 1.5, 0.5)) +
  coord_cartesian(xlim = c(-1, 1.6)) +
  labs(
    x = expression("5-year BMI slope (kg/" * m^2 * "/year)"), y = NULL,
    title = "Recent BMI slope: median, interquartile range, central 95%"
  ) +
  paper_theme

figure <- age_plot / slope_plot +
  plot_layout(heights = c(1, 2.3)) +
  plot_annotation(tag_levels = "A")
save_figure(figure, "design_slope_distributions.png", 7, 6.3)

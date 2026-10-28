# Figure 5: incident outcomes by recent BMI slope and the Diabetes Prevention
# Program (DPP) benchmark.
# A: marginal outcomes; B: outcomes adjusted for age and BMI; C: DPP.
# Inputs: three saved ggplot objects (plot_marginal_outcomes.rds,
#   plot_adjusted_outcomes.rds and plot_dpp_comparison.rds), all stored
#   in the folder data/UKB/processed_data.
# Output figure: plots/outcomes_dpp_summary.png
library(ggplot2)
library(patchwork)

project_dir <- getwd()
plot_data_dir <- file.path(project_dir, "data", "UKB", "processed_data")

marginal_outcomes_plot <- readRDS(
  file.path(plot_data_dir, "plot_marginal_outcomes.rds")
)
# Two BMI ticks keep the narrow facets of panel B legible.
adjusted_outcomes_plot <- readRDS(
  file.path(plot_data_dir, "plot_adjusted_outcomes.rds")
) +
  scale_x_continuous(breaks = c(25, 35))
dpp_plot <- readRDS(
  file.path(plot_data_dir, "plot_dpp_comparison.rds")
)

outcome_panels <- wrap_plots(
  list(marginal_outcomes_plot, adjusted_outcomes_plot),
  nrow = 1, widths = c(7, 5)
)

# Raise the small panel text to a legible size for a ~7-inch printed figure.
paper_theme <- theme(
  text = element_text(size = 8),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 9),
  strip.text = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.tag = element_text(face = "bold", size = 11),
  plot.margin = margin(2, 2, 2, 2)
)

outcomes_summary_plot <- outcome_panels / dpp_plot +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(tag_levels = "A") &
  paper_theme

ggsave(
  file.path(project_dir, "plots", "outcomes_dpp_summary.png"),
  plot = outcomes_summary_plot,
  width = 8, height = 7, units = "in", dpi = 300, bg = "white"
)
cat("Saved plots/outcomes_dpp_summary.png\n")

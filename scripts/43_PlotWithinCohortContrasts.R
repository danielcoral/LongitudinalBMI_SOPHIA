# Figure 2: lower versus upper quartile of recent BMI slope at the same BMI in
# ALSPAC, UK Biobank and the Rotterdam Study (the 17 traits measured in all
# three cohorts). Biomarkers are above DXA measures; differences are in
# standard deviations (SD) of each trait.
# Input: data/cross_cohort/within_cohort_contrasts.csv (script 36).
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

base_plot <- ggplot(
  contrasts, aes(bmi_percentile, estimate_sd, colour = cohort)
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
    y = NULL, colour = NULL
  ) +
  paper_theme +
  theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5))
trait_blocks <- list(
  biomarkers_left = common_traits[1:6],
  biomarkers_right = common_traits[7:11],
  dxa_left = common_traits[12:14],
  dxa_right = common_traits[15:17]
)
blocks <- lapply(trait_blocks, function(traits) {
  base_plot + filter(contrasts, trait %in% traits)
})
blocks$biomarkers_left <- blocks$biomarkers_left +
  labs(tag = "A", x = NULL)
blocks$biomarkers_right <- blocks$biomarkers_right +
  labs(tag = " ", x = NULL)
blocks[1:2] <- lapply(blocks[1:2], function(block) {
  block + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
})
blocks$dxa_left <- blocks$dxa_left + labs(tag = "B")
blocks$dxa_right <- blocks$dxa_right + labs(tag = " ")

# Reserve an empty sixth row to align the upper biomarker columns.
# Measure grobs on a null device to avoid creating an extra Rplots.pdf.
grDevices::pdf(NULL)
block_grobs <- lapply(blocks, function(block) {
  ggplotGrob(block + theme(legend.position = "none"))
})
biomarker_right <- block_grobs$biomarkers_right
last_panel <- biomarker_right$layout$t[
  biomarker_right$layout$name == "panel-1-5"
]
biomarker_right <- gtable::gtable_add_rows(
  biomarker_right,
  grid::unit.c(
    calc_element("panel.spacing.y", blocks$biomarkers_right$theme),
    biomarker_right$heights[last_panel]
  ),
  pos = last_panel
)
block_grobs$biomarkers_right <- biomarker_right
block_widths <- do.call(
  grid::unit.pmax, lapply(block_grobs, function(block) block$widths)
)
display_blocks <- lapply(block_grobs, function(block) {
  block$widths <- block_widths
  wrap_elements(full = block)
})
panels <- wrap_plots(display_blocks, ncol = 2, heights = c(5, 3))
y_label <- grid::textGrob(
  "Upper minus lower quartile of slope (SD units)",
  rot = 90, gp = grid::gpar(fontsize = 9)
)
legend <- gtable::gtable_filter(
  ggplotGrob(blocks$biomarkers_left), "guide-box-bottom"
)
figure <- (
  (wrap_elements(y_label) | panels) + plot_layout(widths = c(0.045, 1))
) / wrap_elements(full = legend) +
  plot_layout(heights = c(1, 0.06))
invisible(grDevices::dev.off())
save_figure(figure, "within_cohort_contrasts.png", 7, 7.2)

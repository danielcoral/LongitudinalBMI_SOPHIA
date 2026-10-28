# Rebuild the supplementary composite figures (primary biomarkers, imaging,
# log liver fat) and trajectory densities from saved results, without refits.
# Outputs: plots/primary_biomarkers.png, imaging.png and
#   log_liver_fat_complete.png and the strict/second-visit density PNGs.
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)
processed_dir <- file.path(getwd(), "data", "UKB", "processed_data")
source(file.path("src", "within_cohort_helpers.R"))
source(file.path("src", "figure_helpers.R"))
source(file.path("src", "trajectory_figure_helpers.R"))
main_theme <- paper_theme

# Raise the small panel text to a legible size for a ~7-inch printed figure.
paper_theme <- theme(
  text = element_text(size = 8),
  axis.text = element_text(size = 8),
  axis.text.x = element_text(size = 6),
  axis.title = element_text(size = 9),
  axis.title.x = element_text(size = 8),
  strip.text = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.tag = element_text(face = "bold", size = 11)
)

primary <- readRDS(file.path(
  processed_dir, "figure_panels_05_GAMRec.rds"
))

# Supplementary Figure 2: slope distribution and variance explained. The
# +0.5 vs -0.5 slope-difference panels are no longer shown, so the
# evaluation-grid markers in panel A (vertical/horizontal lines and points)
# are dropped; only the density regions and the linear trend remain.
density_plot <- primary$bivardistrib_plot
density_plot$layers <- density_plot$layers[1:2]
density_plot <- density_plot +
  facet_wrap(~sex, nrow = 1) +
  scale_x_continuous(breaks = seq(20, 50, 10)) +
  labs(
    x = "BMI at recruitment (kg/m²)",
    y = "5-year BMI slope (kg/m²/year)"
  ) +
  main_theme

# Traits follow the Figure 2 order (liver, lipids, BP, glucose, CRP,
# creatinine), listed top to bottom in panels B-D.
trait_order <- rev(common_traits[1:11])
share_colours <- c(
  Slope = "#D55E00", BMI = "#0072B2", Age = "grey60", Covariates = "grey85"
)
share_scale <- scale_fill_manual(
  values = share_colours, limits = names(share_colours), name = NULL
)
trait_axis <- scale_y_discrete(labels = trait_labels)
bar_axis <- scale_x_continuous(expand = expansion(mult = c(0, 0.04)))

varex_plot <- primary$varex_plot$data |>
  mutate(trait = factor(trait, trait_order)) |>
  ggplot(aes(100 * veval, trait)) +
  geom_col(fill = "grey40", width = 0.7) +
  facet_wrap(~sex, nrow = 1) +
  trait_axis + bar_axis +
  labs(x = "Variance explained (%)", y = NULL) +
  main_theme

# Stacks start with slope at zero so its small share is easy to compare.
stack_order <- rev(names(share_colours))
fai_plot <- primary$fai_plot$data |>
  mutate(
    trait = factor(trait, trait_order),
    term = factor(term, stack_order)
  ) |>
  ggplot(aes(FAI, trait, fill = term)) +
  geom_col(width = 0.7) +
  facet_wrap(~sex, nrow = 1) +
  share_scale + trait_axis + bar_axis +
  labs(x = "Share of explained variance (%)", y = NULL) +
  main_theme

relfai_plot <- primary$relfai_plot$data |>
  mutate(
    trait = factor(trait, trait_order),
    term = factor(term, stack_order)
  ) |>
  ggplot(aes(prop, trait, fill = term)) +
  geom_col(width = 0.7) +
  facet_wrap(~sex, nrow = 1) +
  share_scale + trait_axis + bar_axis +
  labs(x = "BMI vs slope, relative share (%)", y = NULL) +
  guides(fill = "none") +
  main_theme

primary_figure <- (density_plot | varex_plot) /
  ((fai_plot | relfai_plot) + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    panel.spacing.x = grid::unit(1, "lines")
  )
ggsave(
  "plots/primary_biomarkers.png", primary_figure,
  width = 7, height = 6.5, dpi = 300, bg = "white"
)
imaging <- readRDS(file.path(
  processed_dir, "figure_panels_12_GAMImg2i.rds"
))
# Narrow sex facets need fewer ticks as well as smaller x-axis labels.
imaging$varexp_plot <- imaging$varexp_plot +
  scale_x_continuous(breaks = c(0, 0.3, 0.6), labels = \(x) x * 100)
imaging$fai_plot <- imaging$fai_plot +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  labs(y = "Share of explained\nvariance (%)")
imaging$relfai_plot <- imaging$relfai_plot +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  labs(y = "BMI vs slope, relative\nshare (%)")
imaging$pdiff_plot <- imaging$pdiff_plot +
  scale_x_continuous(breaks = c(25, 35))
imaging_figure <- wrap_plots(
  A = imaging$varexp_plot, B = imaging$fai_plot,
  C = imaging$relfai_plot, D = imaging$pdiff_plot,
  design = "ADDD\nBDDD\nCDDD", guides = "keep"
) + plot_annotation(tag_levels = "A") & paper_theme
ggsave(
  "plots/imaging.png", imaging_figure,
  width = 7, height = 5, dpi = 300, bg = "white"
)
log_liver <- c(
  readRDS(file.path(processed_dir, "log_liver_panels.rds")),
  readRDS(file.path(processed_dir, "log_liver_subsampling_panels.rds"))
)
log_liver$C <- log_liver$C +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  labs(y = "Share of explained variance (%)")
for (panel in c("D", "E")) {
  log_liver[[panel]] <- log_liver[[panel]] +
    scale_x_continuous(breaks = c(25, 35)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
log_liver_figure <- wrap_plots(log_liver, ncol = 2) +
  plot_annotation(tag_levels = "A") & paper_theme
ggsave(
  "plots/log_liver_fat_complete.png", log_liver_figure,
  width = 7, height = 5.6, dpi = 300, bg = "white"
)

# Reuse the saved trajectory estimates; do not rerun scripts 09/10's models.
trajectory_figures <- list(
  strict_followup = list(
    input = "lmmtabrec_strict.rds",
    label = "BMI at recruitment (kg/m²)"
  ),
  second_visit = list(
    input = "lmmtab1i.rds",
    label = "BMI at the second visit (kg/m²)"
  )
)
for (name in names(trajectory_figures)) {
  settings <- trajectory_figures[[name]]
  trajectories <- readRDS(file.path(processed_dir, settings$input))
  figure <- trajectory_density_plot(trajectories, settings$label)
  ggsave(
    file.path("plots", paste0(name, "_trajectory_distribution.png")),
    figure, width = 6, height = 3, dpi = 300, bg = "white"
  )
}

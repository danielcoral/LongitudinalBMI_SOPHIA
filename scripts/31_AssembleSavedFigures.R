# Rebuild the supplementary composite figures (primary biomarkers, imaging,
# log liver fat) from saved panels without fitting models.
# Outputs: plots/primary_biomarkers.png, imaging.png and
#   log_liver_fat_complete.png.
library(ggplot2)
library(patchwork)
processed_dir <- file.path(getwd(), "data", "UKB", "processed_data")

# Raise the small panel text to a legible size for a ~7-inch printed figure.
paper_theme <- theme(
  text = element_text(size = 8),
  axis.text = element_text(size = 8),
  axis.title = element_text(size = 9),
  strip.text = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.tag = element_text(face = "bold", size = 11)
)

primary <- readRDS(file.path(
  processed_dir, "figure_panels_05_GAMRec.rds"
))
primary_figure <- wrap_plots(
  A = primary$bivardistrib_plot, B = primary$varex_plot,
  C = primary$fai_plot, D = primary$relfai_plot, E = primary$pdiff_plot,
  design = "ACCDD\nAEEEE\nBEEEE\nBEEEE", guides = "keep"
) + plot_annotation(tag_levels = "A") & paper_theme
ggsave(
  "plots/primary_biomarkers.png", primary_figure,
  width = 7, height = 5, dpi = 300, bg = "white"
)
imaging <- readRDS(file.path(
  processed_dir, "figure_panels_12_GAMImg2i.rds"
))
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
log_liver_figure <- wrap_plots(log_liver, ncol = 2) +
  plot_annotation(tag_levels = "A") & paper_theme
ggsave(
  "plots/log_liver_fat_complete.png", log_liver_figure,
  width = 7, height = 5.6, dpi = 300, bg = "white"
)

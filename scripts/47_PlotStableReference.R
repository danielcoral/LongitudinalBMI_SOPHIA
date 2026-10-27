# Figure 5 and its supplementary versions: gain, stable and loss contrasts in
# UK Biobank and the Rotterdam Study on the common grid (age 65).
# A: gain versus stable (+0.5 vs 0 kg/m2/year); B: stable versus loss (0 vs
# -0.5); C: asymmetry = A minus B. One figure per BMI: 32.5 (main figure),
# 27.5 and 22.5 (supplementary figures).
# Input: data/cross_cohort/stable_contrasts_common_grid.csv (script 40).
# Outputs: plots/stable_reference_bmi32_5.png, stable_reference_bmi27_5.png
#   and stable_reference_bmi22_5.png.
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "figure_helpers.R"))
# Order of Figure 2, with log liver fat after GGT.
trait_order <- c(
  "ALT", "AST", "GGT", "logLiverF_perc", "HDL", "LDL", "TG", "SBP", "DBP",
  "GLU", "CRP", "SCR", "BF_perc", "AndF_perc", "GynF_perc",
  "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
)
stopifnot(setequal(trait_order, names(trait_labels)))
common <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "stable_contrasts_common_grid.csv"),
  show_col_types = FALSE
) |>
  mutate(
    trait = factor(trait, rev(trait_order)),
    cohort = factor(cohort, c("UKB", "RS"))
  )

base_plot <- ggplot(common, aes(estimate_sd, trait, colour = cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(
    aes(xmin = conf_low_sd, xmax = conf_high_sd),
    position = position_dodge(width = 0.7), size = 0.15, linewidth = 0.3
  ) +
  facet_grid(sex ~ .) +
  scale_y_discrete(labels = trait_labels) +
  scale_colour_manual(values = cohort_colors, labels = cohort_labels) +
  labs(y = NULL, colour = NULL) +
  paper_theme
no_strip <- theme(
  strip.text.y = element_blank(), strip.background.y = element_blank()
)
no_y_axis <- theme(
  axis.text.y = element_blank(), axis.ticks.y = element_blank()
)

for (bmi_value in c(32.5, 27.5, 22.5)) {
  shown <- filter(common, bmi == bmi_value)
  gain <- filter(shown, comparison == "gain_vs_stable")
  loss <- filter(shown, comparison == "stable_vs_loss")
  asymmetry <- filter(shown, comparison == "asymmetry")
  # Panels A and B share one x scale so that the gain and loss sides can be
  # compared; the asymmetry has its own scale.
  shared <- range(
    gain$conf_low_sd, gain$conf_high_sd, loss$conf_low_sd, loss$conf_high_sd
  )
  plot_a <- base_plot + gain + no_strip +
    scale_x_continuous(limits = shared) +
    labs(title = "Gain vs stable", x = "Difference (SD units)")
  plot_b <- base_plot + loss + no_y_axis + no_strip +
    scale_x_continuous(limits = shared) +
    labs(title = "Stable vs loss", x = "Difference (SD units)")
  plot_c <- base_plot + asymmetry + no_y_axis +
    labs(title = "Asymmetry", x = "A minus B (SD units)")
  figure <- (plot_a | plot_b | plot_c) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Age 65 years, BMI ", bmi_value, " kg/m\u00b2"),
      tag_levels = "A",
      theme = theme(plot.title = element_text(size = 9, face = "bold"))
    ) &
    theme(legend.position = "bottom")
  save_figure(
    figure, paste0("stable_reference_bmi", gsub("[.]", "_", bmi_value), ".png"),
    7, 8.5
  )
}

# Figure 4: stable-reference contrasts at BMI 30 kg/m2.
# The asymmetry panel is the gain-versus-stable difference minus the
# stable-versus-loss difference.
# Input: data/cross_cohort/stable_contrasts_common_grid.csv (script 39).
# Output figure: stable_reference.png in plots/.
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

gain <- filter(common, comparison == "gain_vs_stable")
loss <- filter(common, comparison == "stable_vs_loss")
asymmetry <- filter(common, comparison == "asymmetry")
# Panels A and B share one x scale so that the gain and loss sides can be
# compared; the asymmetry has its own scale.
shared <- range(
  gain$conf_low_sd, gain$conf_high_sd, loss$conf_low_sd, loss$conf_high_sd
)
plot_a <- base_plot + gain + no_strip +
  scale_x_continuous(limits = shared) +
  labs(title = "a. Gain vs stable", x = "Difference (SD units)")
plot_b <- base_plot + loss + no_y_axis + no_strip +
  scale_x_continuous(limits = shared) +
  labs(title = "b. Stable vs loss", x = "Difference (SD units)")
plot_c <- base_plot + asymmetry + no_y_axis +
  labs(title = "Asymmetry", x = "a minus b (SD units)")
figure <- (plot_a | plot_b | plot_c) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_figure(figure, "stable_reference.png", 7, 8.5)

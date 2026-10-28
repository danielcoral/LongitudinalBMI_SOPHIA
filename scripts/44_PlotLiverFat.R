# Figure 3: regional fat, liver fat and waist-to-hip ratio by recent BMI slope.
# Inputs: saved within-cohort, liver-fat and WHR contrasts.
# Rebuilds the main-text figure from saved cohort, liver-fat and WHR contrasts.
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(patchwork)
library(tidyr)

project_dir <- getwd()
source(file.path(project_dir, "src", "within_cohort_helpers.R"))
source(file.path(project_dir, "src", "figure_helpers.R"))

read_result <- function(path) {
  readr::read_csv(
    file.path(project_dir, "data", "UKB", "results", path),
    show_col_types = FALSE
  )
}

within <- readr::read_csv(
  file.path(project_dir, "data", "cross_cohort",
            "within_cohort_contrasts.csv"),
  show_col_types = FALSE
)
liver <- read_result("liver_fat_quartile_contrasts.csv")
whr <- read_result("whr_trajectory_contrasts.csv")

bmi_levels <- c("p25", "p50", "p75")
cohort_levels <- c("ALSPAC", "UKB", "RS")
cohort_labels <- c("ALSPAC", "UK Biobank", "Rotterdam Study")
cohort_colours <- c(
  ALSPAC = "#E69F00", UKB = "#0072B2", RS = "#009E73"
)

regional <- within |>
  filter(
    in_all_cohorts, trait %in% c("AndF_perc", "GynF_perc")
  ) |>
  select(
    cohort, sex, bmi_percentile, trait,
    estimate_sd, conf_low_sd, conf_high_sd
  ) |>
  pivot_wider(
    names_from = trait,
    values_from = c(estimate_sd, conf_low_sd, conf_high_sd),
    names_glue = "{.value}_{trait}"
  ) |>
  mutate(
    cohort = factor(cohort, cohort_levels, cohort_labels),
    bmi_percentile = factor(bmi_percentile, bmi_levels)
  ) |>
  arrange(cohort, sex, bmi_percentile)

plot_a <- ggplot(
  regional,
  aes(estimate_sd_AndF_perc, estimate_sd_GynF_perc)
) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    colour = "grey50", linewidth = 0.35
  ) +
  geom_path(
    aes(group = interaction(cohort, sex)),
    colour = "grey65", linewidth = 0.35
  ) +
  geom_segment(
    aes(
      x = conf_low_sd_AndF_perc, xend = conf_high_sd_AndF_perc,
      y = estimate_sd_GynF_perc, yend = estimate_sd_GynF_perc
    ),
    colour = "grey55", linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = estimate_sd_AndF_perc, xend = estimate_sd_AndF_perc,
      y = conf_low_sd_GynF_perc, yend = conf_high_sd_GynF_perc
    ),
    colour = "grey55", linewidth = 0.3
  ) +
  geom_point(aes(colour = bmi_percentile), size = 1.8) +
  facet_grid(sex ~ cohort) +
  scale_colour_manual(
    values = c(p25 = "#0072B2", p50 = "#666666", p75 = "#D55E00"),
    labels = c("25th", "50th", "75th"),
    name = "BMI percentile"
  ) +
  coord_equal(xlim = c(-0.25, 0.35), ylim = c(-0.25, 0.35)) +
  labs(
    x = "Android-fat contrast (SD)",
    y = "Gynoid-fat contrast (SD)",
    tag = "A"
  ) +
  paper_theme

liver_age <- liver |>
  filter(grid == "age_bmi_percentile") |>
  mutate(
    sex = factor(sex, c("Female", "Male")),
    bmi_percentile = factor(bmi_percentile, bmi_levels),
    age = factor(age, c(55, 60, 65))
  )

plot_b <- ggplot(
  liver_age,
  aes(
    bmi_percentile, estimate_sd, colour = age,
    group = age
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(
    aes(ymin = conf_low_sd, ymax = conf_high_sd),
    position = position_dodge(width = 0.55),
    linewidth = 0.35
  ) +
  facet_grid(. ~ sex) +
  scale_colour_manual(
    values = c("55" = "#56B4E9", "60" = "#0072B2", "65" = "#D55E00"),
    name = "Age (years)"
  ) +
  scale_x_discrete(labels = c("25th", "50th", "75th")) +
  labs(
    x = "BMI percentile within the imaging sample",
    y = "Liver-fat contrast (SD)",
    tag = "B"
  ) +
  paper_theme

whr_plot <- whr |>
  filter(visit %in% c("rec", "2i")) |>
  mutate(
    visit = factor(
      visit, c("rec", "2i"),
      c("Recruitment", "Imaging visit")
    ),
    sex = factor(sex, c("Female", "Male")),
    bmi_percentile = factor(bmi_percentile, bmi_levels)
  )

plot_c <- ggplot(
  whr_plot,
  aes(
    bmi_percentile, estimate_sd, colour = visit,
    group = visit
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_pointrange(
    aes(ymin = conf_low_sd, ymax = conf_high_sd),
    position = position_dodge(width = 0.55),
    linewidth = 0.35
  ) +
  facet_grid(. ~ sex) +
  scale_colour_manual(
    values = c("Recruitment" = "#0072B2", "Imaging visit" = "#D55E00"),
    name = "WHR assessment"
  ) +
  scale_x_discrete(labels = c("25th", "50th", "75th")) +
  labs(
    x = "BMI percentile within visit and sex",
    y = "WHR contrast (SD)",
    tag = "C"
  ) +
  paper_theme

figure <- (plot_a / plot_b / plot_c) +
  plot_layout(heights = c(1.3, 1, 1)) +
  plot_annotation(tag_levels = NULL)

save_figure(figure, "liver_fat_quartile_contrasts.png", 7, 9.5)

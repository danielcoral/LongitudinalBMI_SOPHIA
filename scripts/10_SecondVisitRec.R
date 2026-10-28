library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(Matrix, warn.conflicts = FALSE)
library(lme4)
library(ggplot2)
library(mgcv)
library(patchwork)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

#-------------------------------------------------------------------------------
# Estimate BMI trajectories up to the second UK Biobank visit (instance 1)

bmilong1i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilong1i.tsv"),
  show_col_types = FALSE
)

second_visit_individual_summary <- bmilong1i |>
  summarise(
    sex = first(sex),
    n_measurements = n(),
    follow_up_years = max(abs(time)),
    .by = eid
  )

second_visit_sample_summary <- second_visit_individual_summary |>
  summarise(
    n_participants = n(),
    n_measurements_median = median(n_measurements),
    n_measurements_q25 = quantile(n_measurements, 0.25),
    n_measurements_q75 = quantile(n_measurements, 0.75),
    proportion_three_or_more = 100 * mean(n_measurements >= 3),
    follow_up_median = median(follow_up_years),
    follow_up_q25 = quantile(follow_up_years, 0.25),
    follow_up_q75 = quantile(follow_up_years, 0.75),
    .by = sex
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

second_visit_overall_summary <- second_visit_individual_summary |>
  summarise(
    n_participants = n(),
    n_measurements_median = median(n_measurements),
    n_measurements_q25 = quantile(n_measurements, 0.25),
    n_measurements_q75 = quantile(n_measurements, 0.75),
    proportion_three_or_more = 100 * mean(n_measurements >= 3),
    follow_up_median = median(follow_up_years),
    follow_up_q25 = quantile(follow_up_years, 0.25),
    follow_up_q75 = quantile(follow_up_years, 0.75)
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

lmmtab1i <- bmilong1i |>
  nest(LongD = -sex) |>
  mutate(
    LongMod = map(
      LongD,
      ~lmer(
        bmi ~ 0 + offset(bmi0) + time + (0 + time | eid),
        data = .x,
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
      )
    ),
    indiv_est = map(LongMod, ~coef(.x)$eid),
    indiv_est = map(indiv_est, data.frame),
    indiv_est = map(indiv_est, setNames, nm = "Slope"),
    indiv_est = map(
      indiv_est,
      ~.x |>
        tibble::rownames_to_column(var = "eid") |>
        mutate(eid = as.numeric(eid))
    )
  )

saveRDS(
  lmmtab1i,
  file.path(projfld, "data", "UKB", "processed_data", "lmmtab1i.rds")
)

average_trajectory_1i <- lmmtab1i |>
  transmute(
    sex,
    model_summary = map(
      LongMod,
      ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE)
    )
  ) |>
  unnest(model_summary) |>
  filter(term == "time") |>
  select(sex, estimate, conf.low, conf.high) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

individual_trajectory_1i <- lmmtab1i |>
  select(sex, indiv_est) |>
  unnest(indiv_est) |>
  summarise(
    slope_median = median(Slope),
    slope_q25 = quantile(Slope, 0.25),
    slope_q75 = quantile(Slope, 0.75),
    .by = sex
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

write_lines(
  c(
    "# Second-visit trajectory sample: overall summary",
    capture.output(print(second_visit_overall_summary, row.names = FALSE)),
    "",
    "# Second-visit trajectory sample: summary by sex",
    capture.output(print(second_visit_sample_summary, row.names = FALSE)),
    "",
    "# Average BMI trajectory up to the second visit (kg/m2/year; 95% CI)",
    capture.output(print(average_trajectory_1i, row.names = FALSE)),
    "",
    "# Individual BMI trajectory distribution by sex (median and IQR)",
    capture.output(print(individual_trajectory_1i, row.names = FALSE))
  ),
  file.path(
    projfld, "data", "UKB", "results",
    "results_second_visit_trajectory.txt"
  )
)
write_csv(
  second_visit_sample_summary,
  file.path(
    projfld, "data", "UKB", "results",
    "table_second_visit_sample_summary.csv"
  )
)
write_csv(
  average_trajectory_1i,
  file.path(
    projfld, "data", "UKB", "results",
    "table_second_visit_average_trajectory.csv"
  )
)

# Supplementary Figure 17: density regions and regression line.
source(file.path(projfld, "src", "trajectory_figure_helpers.R"))
second_visit_trajectory_plot <- trajectory_density_plot(
  lmmtab1i, "BMI at the second visit (kg/m²)"
)

ggsave(
  file.path(
    projfld, "plots",
    "second_visit_trajectory_distribution.png"
  ),
  plot = second_visit_trajectory_plot,
  width = 6, height = 3, units = "in", dpi = 300, bg = "white"
)

#-------------------------------------------------------------------------------
# Biomarker models anchored at the second visit

bm1i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bm1i.tsv"),
  show_col_types = FALSE
)
cv1i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cv1i.tsv"),
  show_col_types = FALSE
)
ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  show_col_types = FALSE
) |>
  mutate(ancestry = factor(ancestry))

covars <- c(
  "CurrSmoke", "Diabetes", "Insulin", "LipidLower", "HTMed",
  "TSI", "CollegeEd", "Menopause", "ancestry"
)
bmrks <- c(
  "ALT", "AST", "GGT", "CRP", "SBP", "DBP",
  "HDL", "LDL", "TG", "GLU", "SCR"
)
smterms <- c(
  "s(age, bs = \"cr\")",
  "s(bmi, bs = \"cr\")",
  "ti(age, bmi)",
  "s(Slope, bs = \"cr\")",
  "ti(bmi, Slope)",
  "ti(Slope, age)"
)

biomarker_data_1i <- bm1i |>
  inner_join(cv1i, by = "eid") |>
  inner_join(ancestrytab, by = "eid") |>
  select(eid, sex, age, bmi, all_of(covars), all_of(bmrks)) |>
  pivot_longer(
    -c(eid, sex, age, bmi, all_of(covars)),
    names_to = "trait",
    values_drop_na = TRUE
  )

slopebmregdf_1i <- biomarker_data_1i |>
  nest(data_by_sex = -sex) |>
  inner_join(lmmtab1i, by = "sex") |>
  transmute(
    sex,
    data_by_sex = map2(indiv_est, data_by_sex, inner_join, by = "eid"),
    data_by_sex = map(data_by_sex, nest, bmdf_1i = -trait)
  ) |>
  unnest(data_by_sex) |>
  mutate(
    mod1_1i = map(
      bmdf_1i,
      ~bam(
        reformulate(
          termlabels = c(covars, smterms),
          response = "value"
        ),
        data = .x,
        drop.unused.levels = TRUE
      )
    )
  )

saveRDS(
  slopebmregdf_1i,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "slopebmregdf_1i.rds"
  )
)

# Compare current second-visit models with the current recruitment models.
slopebmregdf <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf.rds")
)

expectdiffdf <- expand_grid(
  CurrSmoke = 0.08,
  Diabetes = 0.05,
  LipidLower = 0.16,
  HTMed = 0.25,
  Insulin = 0.5,
  TSI = -1.3,
  CollegeEd = 0.3,
  Menopause = 0.4,
  ancestry = "EUR",
  age = c(45, 55, 65),
  bmi = c(22.5, 27.5, 32.5, 37.5),
  Slope = NA
)

second_visit_comparison <- slopebmregdf |>
  select(sex, trait, mod1, bmdf) |>
  inner_join(
    select(slopebmregdf_1i, sex, trait, mod1_1i),
    by = c("sex", "trait")
  ) |>
  transmute(
    sex, trait,
    modpred_recruitment = map(
      mod1,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    modpred_secondvisit = map(
      mod1_1i,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    threshold = map_dbl(bmdf, ~0.1 * sd(.x$value))
  ) |>
  pivot_longer(
    -c(sex, trait, threshold),
    names_to = c(".value", "anchor"),
    names_sep = "_"
  ) |>
  unnest(modpred)

write_csv(
  second_visit_comparison,
  file.path(
    projfld, "data", "UKB", "results",
    "table_second_visit_biomarker_comparison.csv"
  )
)

group_data <- bind_rows(
  tibble(group = "Liver", trait = c("ALT", "AST", "GGT"), panel = 1),
  tibble(group = "Lipids", trait = c("HDL", "LDL", "TG"), panel = 1),
  tibble(group = "BP", trait = c("SBP", "DBP"), panel = 2),
  tibble(group = "Glucose", trait = "GLU", panel = 2),
  tibble(group = "Inflam", trait = "CRP", panel = 2),
  tibble(group = "Kidney", trait = "SCR", panel = 2)
)

second_visit_biomarker_plot <- second_visit_comparison |>
  inner_join(group_data, by = "trait") |>
  mutate(
    group = factor(
      group,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney", "Adipose"
      )
    ),
    anchor = factor(
      anchor,
      levels = c("recruitment", "secondvisit"),
      labels = c("Recruitment", "Second visit")
    )
  ) |>
  nest(plot_data = -panel) |>
  mutate(
    plot = map(
      plot_data,
      ~.x |>
        ggplot(aes(bmi, estimate)) +
        geom_rect(
          aes(ymin = -threshold, ymax = threshold),
          xmin = -Inf, xmax = Inf,
          fill = "grey70", alpha = 0.3
        ) +
        geom_hline(yintercept = 0, linetype = "11", linewidth = 0.2) +
        geom_hline(
          aes(yintercept = -threshold),
          linetype = "11", linewidth = 0.2
        ) +
        geom_hline(
          aes(yintercept = threshold),
          linetype = "11", linewidth = 0.2
        ) +
        geom_linerange(
          aes(ymin = conf.low, ymax = conf.high, group = anchor),
          position = position_dodge(width = 1.25),
          linewidth = 0.3
        ) +
        geom_point(
          aes(group = anchor, color = anchor),
          position = position_dodge(width = 1.25),
          size = 0.5
        ) +
        ggh4x::facet_nested(
          group + trait ~ sex + age,
          scales = "free",
          strip = ggh4x::strip_nested(size = "variable")
        ) +
        labs(
          x = "BMI",
          y = "RG - RL difference",
          color = "Anchor visit"
        ) +
        theme_bw()
    ),
    plot = map2(
      panel, plot,
      ~if (.x == 2) .y + theme(axis.title.y = element_blank()) else .y
    )
  ) |>
  pull(plot) |>
  wrap_plots(nrow = 1, guides = "collect") &
  theme(
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = -10, b = -10),
    legend.key.width = unit(1, "mm")
  )

ggsave(
  file.path(
    projfld, "plots",
    "second_visit_biomarker_comparison.png"
  ),
  plot = second_visit_biomarker_comparison_plot,
  width = 7, height = 4, units = "in", dpi = 300, bg = "white"
)

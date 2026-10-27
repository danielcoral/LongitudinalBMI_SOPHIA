library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(Matrix, warn.conflicts = FALSE)
library(lme4)
library(ggplot2)
library(mgcv)
library(patchwork)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

#-------------------------------------------------------------------------------
# Select participants with at least three well-spaced BMI measurements

bmilongrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv"),
  show_col_types = FALSE
)

strict_ids <- bmilongrec |>
  summarise(
    n_measurements = n(),
    has_baseline = any(time == 0),
    has_mid_measurement = any(time >= -2.5 & time <= -0.5),
    has_early_measurement = any(time < -2.5),
    .by = eid
  ) |>
  filter(
    n_measurements >= 3,
    has_baseline,
    has_mid_measurement,
    has_early_measurement
  ) |>
  select(eid)

bmilongrec_strict <- bmilongrec |>
  semi_join(strict_ids, by = "eid")

strict_sample_summary <- bmilongrec_strict |>
  summarise(
    sex = first(sex),
    n_measurements = n(),
    follow_up_years = max(abs(time)),
    .by = eid
  ) |>
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

write_tsv(
  bmilongrec_strict,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "bmilongrec_strict.tsv"
  )
)

#-------------------------------------------------------------------------------
# Estimate trajectories in the strict-follow-up sample

lmmtabrec_strict <- bmilongrec_strict |>
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
  lmmtabrec_strict,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "lmmtabrec_strict.rds"
  )
)

average_trajectory_strict <- lmmtabrec_strict |>
  transmute(
    sex,
    estimate = map(
      LongMod,
      ~broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE)
    )
  ) |>
  unnest(estimate) |>
  filter(term == "time") |>
  select(sex, estimate, conf.low, conf.high) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

write_lines(
  c(
    "# Strict-follow-up criteria",
    "At least 3 BMI measurements including baseline (time 0),",
    "at least one measurement 0.5 to 2.5 years before recruitment,",
    "and at least one measurement more than 2.5 years before recruitment.",
    "",
    "# Strict-follow-up sample summary",
    capture.output(print(strict_sample_summary, row.names = FALSE)),
    "",
    "# Average BMI trajectory (kg/m2/year; estimate and 95% CI)",
    capture.output(print(average_trajectory_strict, row.names = FALSE))
  ),
  file.path(
    projfld, "data", "UKB", "results",
    "results_strict_followup_trajectory.txt"
  )
)
write_csv(
  strict_sample_summary,
  file.path(
    projfld, "data", "UKB", "results",
    "table_strict_followup_sample_summary.csv"
  )
)
write_csv(
  average_trajectory_strict,
  file.path(
    projfld, "data", "UKB", "results",
    "table_strict_followup_average_trajectory.csv"
  )
)

# Supplementary Figure 9
strict_trajectory_plot <- lmmtabrec_strict |>
  transmute(
    sex,
    plot_data = map(LongD, select, eid, bmi0),
    plot_data = map(plot_data, distinct),
    plot_data = map2(plot_data, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(plot_data) |>
  ggplot(aes(x = bmi0, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
  geom_vline(
    xintercept = c(22.5, 27.5, 32.5, 37.5, 42.5),
    linetype = "11", color = "magenta"
  ) +
  geom_hline(
    yintercept = c(-0.5, 0, 0.5),
    linetype = "11", color = "darkgreen"
  ) +
  geom_point(
    data = expand_grid(
      bmi0 = c(22.5, 27.5, 32.5, 37.5),
      Slope = c(-0.5, 0, 0.5)
    ),
    color = "blue"
  ) +
  facet_wrap(~sex) +
  theme_bw() +
  labs(
    x = "BMI at recruitment (kg/m²)",
    y = "Individual slope (kg/m² per year)",
    alpha = "Probability\ndensity"
  ) +
  theme(
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  coord_cartesian(xlim = c(20, 45), ylim = c(-2.5, 2.5))

ggsave(
  file.path(
    projfld, "plots",
    "strict_followup_trajectory_distribution.png"
  ),
  plot = strict_trajectory_plot,
  width = 6, height = 3, units = "in", dpi = 300, bg = "white"
)

#-------------------------------------------------------------------------------
# Biomarker models in the strict-follow-up sample

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)
cvrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv"),
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

biomarker_data <- bmrec |>
  inner_join(cvrec, by = "eid") |>
  inner_join(ancestrytab, by = "eid") |>
  select(eid, sex, age, bmi, all_of(covars), all_of(bmrks)) |>
  pivot_longer(
    -c(eid, sex, age, bmi, all_of(covars)),
    names_to = "trait",
    values_drop_na = TRUE
  )

slopebmregdf_strict <- biomarker_data |>
  nest(data_by_sex = -sex) |>
  inner_join(lmmtabrec_strict, by = "sex") |>
  transmute(
    sex,
    data_by_sex = map2(indiv_est, data_by_sex, inner_join, by = "eid"),
    data_by_sex = map(data_by_sex, nest, bmdf = -trait)
  ) |>
  unnest(data_by_sex) |>
  mutate(
    mod_strict = map(
      bmdf,
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
  slopebmregdf_strict,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "slopebmregdf_strict.rds"
  )
)

# Recalculate the original estimates from the current primary models so that the
# comparison always uses the same prediction grid and current pipeline results.
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

strict_comparison <- slopebmregdf_strict |>
  select(sex, trait, mod_strict) |>
  inner_join(
    select(slopebmregdf, sex, trait, mod1, bmdf),
    by = c("sex", "trait")
  ) |>
  transmute(
    sex, trait,
    modpred_all = map(
      mod1,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    modpred_strict = map(
      mod_strict,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    threshold = map_dbl(bmdf, ~0.1 * sd(.x$value))
  ) |>
  pivot_longer(
    -c(sex, trait, threshold),
    names_to = c(".value", "sample"),
    names_sep = "_"
  ) |>
  unnest(modpred)

write_csv(
  strict_comparison,
  file.path(
    projfld, "data", "UKB", "results",
    "table_strict_followup_biomarker_comparison.csv"
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

strict_biomarker_comparison_plot <- strict_comparison |>
  inner_join(group_data, by = "trait") |>
  mutate(
    group = factor(
      group,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney", "Adipose"
      )
    ),
    sample = factor(
      sample,
      levels = c("all", "strict"),
      labels = c("All", "Strict follow-up")
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
          aes(ymin = conf.low, ymax = conf.high, group = sample),
          position = position_dodge(width = 1.25),
          linewidth = 0.3
        ) +
        geom_point(
          aes(group = sample, color = sample),
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
          color = "Sample"
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
    "strict_followup_biomarker_comparison.png"
  ),
  plot = strict_biomarker_comparison_plot,
  width = 7, height = 4, units = "in", dpi = 300, bg = "white"
)

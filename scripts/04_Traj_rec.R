library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(Matrix, warn.conflicts = FALSE)
library(lme4)
library(ggplot2)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

bmilongrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv"),
  show_col_types = FALSE
)

lmmtabrec <- nest(bmilongrec, LongD = -sex)

lmmtabrec <- lmmtabrec |>
  mutate(
    LongMod = map(
      LongD,
      ~lmer(
        formula = bmi ~ 0 + offset(bmi0) + time + (0 + time | eid),
        data = .x,
        ## Additional options for faster output
        REML = FALSE,
        control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
      )
    )
  )

lmmtabrec <- lmmtabrec |>
  mutate(
    ## Extracting estimates
    indiv_est = map(LongMod, ~coef(.x)$eid),
    ## Putting them in a table
    indiv_est = map(indiv_est, data.frame),
    ## Formatting table
    indiv_est = map(indiv_est, setNames, nm = c("Slope")),
    indiv_est = map(indiv_est, tibble::rownames_to_column, var = "eid"),
    indiv_est = map(indiv_est, mutate, eid = as.numeric(eid))
  )

saveRDS(
  lmmtabrec,
  file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
)

average_trajectory <- lmmtabrec |>
  transmute(
    sex,
    LongMod = map(
      LongMod, broom.mixed::tidy,
      effects = "fixed", conf.int = TRUE
    )
  ) |>
  unnest(LongMod) |>
  filter(term == "time") |>
  select(sex, estimate, conf.low, conf.high) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))

individual_trajectory <- lmmtabrec |>
  select(sex, indiv_est) |>
  unnest(indiv_est) |>
  group_by(sex) |>
  summarise(
    across(
      Slope,
      list(
        Med = \(x) median(x, na.rm = TRUE),
        Q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        Q75 = \(x) quantile(x, probs = .75, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)

bv_age <- lmmtabrec |>
  transmute(
    cordat = map(indiv_est, inner_join, bmrec, by = "eid")
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = age, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
  facet_wrap(~sex) +
  theme_bw() +
  labs(
    x = "Age at recruitment (years)",
    y = "Individual slope (kg/m² per year)",
    alpha = "Probability\ndensity"
  ) +
  theme(
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  coord_cartesian(xlim = c(40, 70), ylim = c(-2.5, 2.5))

bv_nm <- lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, count, eid),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid"),
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = log(n), y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
  facet_wrap(~sex) +
  theme_bw() +
  labs(
    x = "Log number of measures",
    y = "Individual slope (kg/m² per year)",
    alpha = "Probability\ndensity"
  ) +
  theme(
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  coord_cartesian(xlim = c(NA, 4), ylim = c(-2.5, 2.5))

bv_fut <- lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, group_by, eid),
    cordat = map(cordat, summarise, maxT = max(abs(time))),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid"),
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = maxT, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
  facet_wrap(~sex) +
  theme_bw() +
  labs(
    x = "Follow-up time",
    y = "Individual slope (kg/m² per year)",
    alpha = "Probability\ndensity"
  ) +
  theme(
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  coord_cartesian(ylim = c(-2.5, 2.5))

trajectory_correlations <- lmmtabrec |>
  transmute(
    sex,
    data = map2(
      LongD, indiv_est,
      \(long_data, slopes) {
        long_data |>
          left_join(select(bmrec, eid, age), by = "eid") |>
          summarise(
            bmi0 = first(bmi0),
            age = first(age),
            n_measurements = n(),
            follow_up_years = max(abs(time)),
            .by = eid
          ) |>
          inner_join(slopes, by = "eid")
      }
    )
  ) |>
  unnest(data) |>
  summarise(
    rho_current_bmi = cor(bmi0, Slope, method = "spearman"),
    rho_age = cor(age, Slope, method = "spearman", use = "complete.obs"),
    rho_log_measurements = cor(log(n_measurements), Slope, method = "spearman"),
    rho_follow_up = cor(follow_up_years, Slope, method = "spearman"),
    .by = sex
  ) |>
  mutate(across(starts_with("rho_"), \(x) round(x, 3)))

trajectory_text <- c(
  "# Average BMI trajectory by sex (kg/m2/year; estimate and 95% CI)",
  capture.output(print(average_trajectory, row.names = FALSE)),
  "",
  "# Individual BMI trajectory distribution by sex (median and IQR)",
  capture.output(print(individual_trajectory, row.names = FALSE)),
  "",
  "# Spearman correlations with individual BMI trajectory",
  capture.output(print(trajectory_correlations, row.names = FALSE))
)
write_lines(
  trajectory_text,
  file.path(projfld, "data", "UKB", "results", "results_primary_trajectory.txt")
)
write_csv(
  average_trajectory,
  file.path(
    projfld, "data", "UKB", "results",
    "table_primary_average_trajectory.csv"
  )
)
write_csv(
  trajectory_correlations,
  file.path(
    projfld, "data", "UKB", "results",
    "table_primary_trajectory_correlations.csv"
  )
)

trajectory_correlates_plot <- patchwork::wrap_plots(
  bv_age, bv_nm, bv_fut, ncol = 1
) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(projfld, "plots", "trajectory_correlates.png"),
  plot = trajectory_correlates_plot,
  width = 6, height = 9, units = "in", dpi = 300,
  bg = "white"
)

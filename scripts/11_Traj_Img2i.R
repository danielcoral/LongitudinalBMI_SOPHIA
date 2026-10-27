library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(Matrix, warn.conflicts = FALSE)
library(lme4)
library(ggplot2)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

bmilong2i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilong2i.tsv"),
  show_col_types = FALSE
)

trajectory_data_summary_2i <- bmilong2i |>
  group_by(eid) |>
  summarise(
    NMs = n(),
    FUT = min(time) - max(time),
    .groups = "drop"
  ) |>
  summarise(
    across(
      c(NMs, FUT),
      list(
        Med = \(x) median(x, na.rm = TRUE),
        q25 = \(x) quantile(x, probs = .25, na.rm = TRUE),
        q75 = \(x) quantile(x, probs = .75, na.rm = TRUE)
      )
    ),
    Prop3Ms = 100 * mean(NMs >= 3, na.rm = TRUE)
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

lmmtab2i <- nest(bmilong2i, LongD = -sex)

lmmtab2i <- lmmtab2i |>
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

lmmtab2i <- lmmtab2i |>
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
  lmmtab2i,
  file.path(projfld, "data", "UKB", "processed_data", "lmmtab2i.rds")
)

average_trajectory_2i <- lmmtab2i |>
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

write_lines(
  c(
    "# Trajectory data before imaging visit",
    capture.output(print(trajectory_data_summary_2i, row.names = FALSE)),
    "",
    "# Average BMI trajectory before imaging by sex (kg/m2/year; 95% CI)",
    capture.output(print(average_trajectory_2i, row.names = FALSE))
  ),
  file.path(projfld, "data", "UKB", "results", "results_imaging_trajectory.txt")
)
write_csv(
  average_trajectory_2i,
  file.path(
    projfld, "data", "UKB", "results",
    "table_imaging_average_trajectory.csv"
  )
)

imaging_trajectory_plot <- lmmtab2i |>
  transmute(
    sex,
    cordat = map(LongD, select, eid, bmi0),
    cordat = map(cordat, unique),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = bmi0, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_vline(
    xintercept = c(22.5, 27.5, 32.5, 37.5, 42.5),
    linetype = "11", color = "magenta"
  ) +
  geom_hline(yintercept = c(-.5, 0, .5), linetype = "11", color = "darkgreen") +
  geom_point(
    data = expand_grid(
      bmi0 = c(22.5, 27.5, 32.5, 37.5),
      Slope = c(-.5, 0, .5)
    ),
    color = "blue"
  ) +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
  facet_wrap(~sex) +
  theme_bw() +
  labs(
    x = "BMI at the imaging visit (kg/m²)",
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
  file.path(projfld, "plots", "imaging_trajectory_distribution.png"),
  plot = imaging_trajectory_plot,
  width = 6, height = 3, units = "in", dpi = 300,
  bg = "white"
)
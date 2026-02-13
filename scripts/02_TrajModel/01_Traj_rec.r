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

head(bmilongrec)

lmmtabrec <- nest(bmilongrec, LongD = -sex)

print(lmmtabrec)

lmmtabrec <- lmmtabrec %>%
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

print(lmmtabrec)

lmmtabrec <- lmmtabrec %>%
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

print(lmmtabrec)

saveRDS(
  lmmtabrec,
  file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
)

lmmtabrec <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
)

lmmtabrec |>
  transmute(
    sex,
    LongMod = map(
      LongMod, broom.mixed::tidy,
      effects = "fixed", conf.int = TRUE
    )
  ) |>
  unnest(LongMod) |>
  select(sex, estimate, conf.low, conf.high) |>
  mutate(across(where(is.numeric), ~round(.x, 2)))

lmmtabrec |>
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
    )
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, digits = 2)))

lmmtabrec |>
  transmute(
    sex,
    cordat = map(indiv_est, select, eid, Slope),
    cordat = map2(cordat, LongD, left_join, by = "eid"),
    cordat = round(
      map_dbl(cordat, with, cor(bmi0, Slope, method = "spearman")),
      2
    )
  ) |>
  unnest(cordat)

options(repr.plot.width = 6, repr.plot.height = 3)
lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, select, eid, bmi0),
    cordat = map(cordat, unique),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = bmi0, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
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

options(repr.plot.width = 6, repr.plot.height = 3)
lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, select, eid, bmi0),
    cordat = map(cordat, unique),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = bmi0, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red") +
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

lmmtabrec |>
  transmute(
    sex,
    LongD = map2(LongD, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(LongD) |>
  filter(bmi0 >= 30, bmi0 <= 40, Slope <= -4) |>
  arrange(eid, time) |>
  add_count(eid) |>
  arrange(n) |>
  ggplot(aes(x = time, y = bmi, group = eid)) +
  geom_line(alpha = 0.3) +
  geom_abline(
    aes(intercept = bmi0, slope = Slope),
    color = "red", alpha = 0.3
  ) +
  facet_wrap(~paste(eid, round(Slope, 2), sep = " - ")) +
  theme_bw() +
  labs(x = "Time since recruitment (years)", y = "BMI (kg/m²)")

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)

head(bmrec)

lmmtabrec |>
  transmute(
    cordat = map(indiv_est, inner_join, bmrec, by = "eid")
  ) |>
  unnest(cordat) |>
  summarise(rho = round(cor(age, Slope, method = "spearman"), 2))

lmmtabrec |>
  transmute(
    sex,
    cordat = map(indiv_est, inner_join, bmrec, by = "eid"),
    cordat = round(
      map_dbl(
        cordat,
        with, cor(age, Slope, method = "spearman")
      ),
      2
    )
  ) |>
  unnest(cordat)

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
bv_age

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, count, eid)
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = n)) +
  geom_histogram(bins = 100) +
  facet_wrap(~sex) +
  theme_bw() +
  labs(x = "Number of measures", y = "N") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  )

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, count, eid),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  summarise(rho = round(cor(Slope, log(n), method = "spearman"), 2))

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, count, eid),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid"),
    cordat = round(
      map_dbl(cordat, with, cor(Slope, log(n), method = "spearman")),
      2
    )
  ) |>
  unnest(cordat)

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
bv_nm

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, group_by, eid),
    cordat = map(cordat, summarise, maxT = max(abs(time)))
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = maxT)) +
  geom_histogram(bins = 100) +
  facet_wrap(~sex) +
  theme_bw() +
  labs(x = "Follow-up time", y = "N") +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  )

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, group_by, eid),
    cordat = map(cordat, summarise, maxT = max(abs(time))),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  summarise(rho = round(cor(maxT, Slope, method = "spearman"), 2))

lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, group_by, eid),
    cordat = map(cordat, summarise, maxT = max(abs(time))),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid"),
    cordat = round(
      map_dbl(cordat, with, cor(Slope, maxT, method = "spearman")),
      2
    )
  ) |>
  unnest(cordat)

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
bv_fut

options(repr.plot.width = 6, repr.plot.height = 9)
patchwork::wrap_plots(
  bv_age, bv_nm, bv_fut, ncol = 1
) +
  patchwork::plot_annotation(tag_levels = "A")

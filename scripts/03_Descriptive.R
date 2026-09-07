# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)

# Project folder
projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

# Longitudinal BMI data until recruitment
bmilongrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv"),
  show_col_types = FALSE
)

# Inclusion indicator
ind_df <- bmilongrec |>
  select(eid) |>
  distinct() |>
  mutate(included = 1)

# Biomarker data at recruitment
bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)

# Covariate data at recruitment
cvrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv"),
  show_col_types = FALSE
)
ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  show_col_types = FALSE
)


# Summary of coviarate data at recruitment by sex
sumcvrec <- cvrec |>

  inner_join(bmrec[, c("eid", "sex")], by = "eid") |>
  summarise(
    across(-eid, \(x) mean(x, na.rm = TRUE)),
    .by = sex
  )
sumcvrec

# Survival data after recruitment
outdat <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv"),
  show_col_types = FALSE
)

# Summary of survival data
sumoutdat <- outdat |>
  inner_join(bmrec[, c("eid", "sex")], by = "eid") |>
  nest(Dat = -c(sex, outcome)) |>
  mutate(
    sex, outcome,
    km_est = map(
      Dat,
      ~{
        sums <- survfit(Surv(time, event) ~ 1, data = .x) |>
          summary(times = 10)
        tibble(
          estimate = sums$surv,
          conf.low = sums$lower,
          conf.high = sums$upper
        ) |>
          mutate(
            across(everything(), ~ 100 * (1 - .x)),
            n_events = sum(.x$event),
            n_py = sum(.x$time),
            inc_rate = 100000 * (n_events / n_py)
          )
      }
    )
  ) |>
  unnest(km_est)

sumoutdat

write_tsv(
  sumoutdat,
  file.path(projfld, "data", "UKB", "results", "sumoutdat.tsv")
)
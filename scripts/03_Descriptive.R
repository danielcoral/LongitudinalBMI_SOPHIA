# Libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)

# Project folder
projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

# Total number of individuals available
totaln <- system(
  paste(
    "wc -l",
    file.path(
      projfld, "data", "UKB", "raw_data",
      "Weight_Recruitment_participant.csv"
    )
  ),
  intern = TRUE
) |>
  stringr::str_extract("^[0-9]+")
descriptext <- paste(
  "Total UKB participants available: ",
  as.numeric(totaln) - 1
)

totalwgp <- system(
  paste(
    "cut -f 1",
    file.path(
      projfld, "data", "UKB", "raw_data",
      "OnlyWeights_gp.tsv"
    ), "| tail -n +2 |",
    "sort -u | wc -l"
  ),
  intern = TRUE
)
descriptext <- c(
  descriptext, "",
  paste("N with GP weight data: ", totalwgp)
)

# Longitudinal BMI data until recruitment
bmilongrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmilongrec.tsv"),
  show_col_types = FALSE
)

# Total number of participants with longitudinal BMI data until recruitment
descriptext <- c(
  descriptext, "",
  "Total N with valid longitudinal BMI data until recruitment:",
  length(unique(bmilongrec$eid))
)

# By sex
nbysex <- bmilongrec |>
  group_by(sex) |>
  summarise(N = n_distinct(eid))

descriptext <- c(
  descriptext, "",
  "N with longitudinal BMI data until recruitment by sex:",
  capture.output(data.frame(nbysex))
)

# Inclusion indicator
ind_df <- bmilongrec |>
  select(eid) |>
  distinct() |>
  mutate(inclstatus = "Included")

nrow(ind_df)

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

# Helper function to expand categorical variables
expand_catvars <- function(df) {
  dfcat <- df |>
    select(
      where(
        ~length(unique(na.omit(.x))) > 2 &&
          length(unique(na.omit(.x))) <= 10
      )
    )
  df <- df |>
    select(-all_of(names(dfcat)))
  expdfcat <- matrix(NA, nrow = nrow(df))
  for (varnm in names(dfcat)) {
    fv <- factor(dfcat[[varnm]])
    mm <- model.matrix.lm(~ fv - 1, na.action = na.pass)
    colnames(mm) <- paste0(varnm, "_", levels(fv))
    expdfcat <- cbind(expdfcat, mm)
  }
  bind_cols(df, data.frame(expdfcat))
}

descriptable <- ind_df |>
  full_join(select(bmrec, -WHR), by = "eid") |>
  full_join(cvrec, by = "eid") |>
  left_join(ancestrytab, by = "eid") |>
  replace_na(list(inclstatus = "Excluded")) |>
  expand_catvars() |>
  pivot_longer(
    cols = -c(eid, sex, inclstatus),
    names_to = "var_nm",
    values_to = "var_val"
  ) |>
  summarise(
    Navail = sum(!is.na(var_val)),
    var_valsum = if_else(
      length(unique(na.omit(var_val))) <= 10,
      paste(round(100 * mean(var_val, na.rm = TRUE), 2), "%"),
      paste0(
        round(mean(var_val, na.rm = TRUE), 2),
        " (",
        round(quantile(var_val, probs = 0.25, na.rm = TRUE), 2),
        " - ",
        round(quantile(var_val, probs = 0.75, na.rm = TRUE), 2),
        ")"
      )
    ),
    .by = c(var_nm, sex, inclstatus)
  ) |>
  pivot_wider(
    names_from = inclstatus,
    values_from = c(Navail, var_valsum),
    names_vary = "slowest"
  )
write_tsv(
  descriptable,
  file.path(projfld, "data", "UKB", "results", "descriptable.tsv")
)

# Follow-up time summary
overall_fut <- bmilongrec |>
  summarise(
    min_time = min(time),
    .by = eid
  ) |>
  summarise(
    across(
      min_time,
      list(
        q25 = ~round(quantile(.x, probs = 0.25), 2),
        median = ~round(median(.x), 2),
        q75 = ~round(quantile(.x, probs = 0.75), 2)
      )
    )
  )
descriptext <- c(
  descriptext, "",
  "Overall follow-up time (years):",
  capture.output(data.frame(overall_fut))
)

fut_bysex <- bmilongrec |>
  summarise(
    sex = first(sex),
    min_time = min(time),
    .by = eid
  ) |>
  summarise(
    across(
      min_time,
      list(
        q25 = ~round(quantile(.x, probs = 0.25), 2),
        median = ~round(median(.x), 2),
        q75 = ~round(quantile(.x, probs = 0.75), 2)
      )
    ),
    .by = sex
  )
descriptext <- c(
  descriptext, "",
  "Follow-up time by sex (years):",
  capture.output(data.frame(fut_bysex))
)

# N measures summary
overall_nmeas <- bmilongrec |>
  count(eid, sex, name = "nmeas") |>
  summarise(
    across(
      nmeas,
      list(
        q25 = ~round(quantile(.x, probs = 0.25), 2),
        median = ~round(median(.x), 2),
        q75 = ~round(quantile(.x, probs = 0.75), 2),
        prop3ormore = ~round(100 * mean(.x >= 3), 2)
      )
    )
  )
descriptext <- c(
  descriptext, "",
  "Overall number of measurements:",
  capture.output(data.frame(overall_nmeas))
)

nmeas_bysex <- bmilongrec |>
  count(eid, sex, name = "nmeas") |>
  summarise(
    across(
      nmeas,
      list(
        q25 = ~round(quantile(.x, probs = 0.25), 2),
        median = ~round(median(.x), 2),
        q75 = ~round(quantile(.x, probs = 0.75), 2)
      )
    ),
    .by = sex
  )
descriptext <- c(
  descriptext, "",
  "Number of measurements by sex:",
  capture.output(data.frame(nmeas_bysex))
)

write_lines(
  descriptext,
  file.path(projfld, "data", "UKB", "results", "descriptext.txt")
)

# Survival data after recruitment
outdat <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv"),
  show_col_types = FALSE
)

# Summary of survival data
sumoutdat <- outdat |>
  semi_join(distinct(bmilongrec, eid), by = "eid") |>
  inner_join(bmrec[, c("eid", "sex")], by = "eid") |>
  summarise(
    n_events = sum(event),
    n_py = sum(time),
    inc_rate = round(100000 * (n_events / n_py), 2),
    .by = c("sex", "outcome")
  )
write_tsv(
  sumoutdat,
  file.path(projfld, "data", "UKB", "results", "sumoutdat.tsv")
)

write_csv(
  descriptable,
  file.path(
    projfld, "data", "UKB", "results",
    "participant_characteristics.csv"
  )
)

trajectory_followup <- bind_rows(
  mutate(nmeas_bysex, variable = "Number of measurements"),
  mutate(
    fut_bysex,
    variable = "Follow-up time (years)",
    nmeas_q25 = min_time_q25,
    nmeas_median = min_time_median,
    nmeas_q75 = min_time_q75
  ) |>
    select(sex, starts_with("nmeas_"), variable)
) |>
  rename(q25 = nmeas_q25, median = nmeas_median, q75 = nmeas_q75)
write_csv(
  trajectory_followup,
  file.path(
    projfld, "data", "UKB", "results",
    "trajectory_followup.csv"
  )
)
write_csv(
  sumoutdat,
  file.path(
    projfld, "data", "UKB", "results",
    "outcome_summary.csv"
  )
)
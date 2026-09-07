# Data wrangling
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(broom.mixed)
library(mgcv)
library(metagam)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA/"

input_table_bysex <- bind_rows(
  "rec" = readRDS(
    file = file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
  ),
  "2i" = readRDS(
    file = file.path(projfld, "data", "UKB", "processed_data", "lmmtab2i.rds")
  ),
  .id = "reftime"
)
print(input_table_bysex)

bmdfs <- tibble(
  reftime = c("rec", "2i"),
) |>
  mutate(
    bmdat = map(
      reftime,
      ~switch(
        .x,
        "rec" = read_tsv(
          file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
          show_col_types = FALSE
        ),
        "2i" = read_tsv(
          file.path(projfld, "data", "UKB", "processed_data", "img2i.tsv"),
          show_col_types = FALSE
        )
      )
    )
  )
print(bmdfs)

covdfs <- tibble(
  reftime = c("rec", "2i"),
) |>
  mutate(
    covardat = map(
      reftime,
      ~switch(
        .x,
        "rec" = read_tsv(
          file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv"),
          show_col_types = FALSE
        ),
        "2i" = read_tsv(
          file.path(projfld, "data", "UKB", "processed_data", "cv2i.tsv"),
          show_col_types = FALSE
        )
      )
    )
  )
print(covdfs)

agebmi <- c("age", "bmi")

bmrks <- c(
  # Biochemistry markers
  "ALT", "AST", "CRP", "DBP", "GGT", "HDL",
  "LDL", "GLU", "SBP", "SCR", "TG", "WHR",
  # DEXA markers
  "BF_perc", "AndF_perc", "GynF_perc",
  "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
)

covars <- c(
  # Current smoking status
  "CurrSmoke",
  # Comorbidities and medications
  "Diabetes", "LipidLower",
  "HTMed", "Insulin",
  # Socioeconomic covariates
  "TSI", "CollegeEd",
  # Sex-specific covariates
  "Menopause"
)

joindf <- input_table_bysex |>
  inner_join(bmdfs, by = "reftime") |>
  inner_join(covdfs, by = "reftime") |>
  mutate(
    modx = map2(indiv_est, bmdat, inner_join, by = "eid"),
    modx = map2(modx, covardat, inner_join, by = "eid"),
    modx = map(
      modx,
      select,
      eid,
      all_of(agebmi),
      Slope,
      any_of(covars),
      any_of(bmrks)
    )
  )
print(joindf)

regresdat <- tibble(
  reftime = c("rec", "2i"),
) |>
  mutate(
    modlist = map(
      reftime,
      ~switch(
        .x,
        "rec" = readRDS(
          file = file.path(
            projfld, "data", "UKB", "processed_data", "slopebmregdf.rds"
          )
        ),
        "2i" = readRDS(
          file = file.path(
            projfld, "data", "UKB", "processed_data", "slopebmregdf_2i.rds"
          )
        )
      )
    ),
    modlist = map(modlist, nest, modres = -sex)
  ) |>
  unnest(modlist) |>
  mutate(
    modres = map(
      modres,
      transmute,
      trait,
      mod1 = map(mod1, strip_rawdata)
    )
  )
print(regresdat)

joindf <- joindf |>
  inner_join(regresdat, by = c("reftime", "sex")) |>
  rename(LongData = LongD)
print(joindf)

joindf <- joindf %>%
  mutate(
    LongDatSum = map(
      LongData,
      ~{
        .x %>%
          # Extract follow-up time and measures
          group_by(eid) %>%
          summarise(T = min(time), NM = n()) %>%
          # Calculate summary statistics
          summarise(
            TotalInd = n(), TotalT = sum(abs(T)), TotalM = sum(NM),
            across(
              c(T, NM),
              list(
                Mn = mean, Sd = sd, Md = median,
                Q25 = \(x) quantile(x, probs = .25),
                Q75 = \(x) quantile(x, probs = .75)
              )
            ),
            NM_mode = tibble(X = NM) %>%
              count(X) %>%
              slice_max(n) %>%
              slice_min(as.numeric(X)) %>%
              pull(X),
            NM_modeprop = 100 * sum(NM == NM_mode) / TotalInd
          )
      }
    )
  )

joindf <- joindf %>%
  mutate(
    # Average change over time
    LongModSum = map(
      LongMod,
      broom.mixed:::tidy.merMod,
      effects = "fixed",
      conf.int = TRUE
    ),
    # Cleaning
    LongModSum = map(LongModSum, filter, term == "time"),
    LongModSum = map(LongModSum, select, -c(effect, term))
  )

joindf <- joindf %>%
  mutate(TotalN = map_dbl(modx, nrow))

# Continous variables, change accordingly
contvars <- c(agebmi, bmrks, "Slope", "TSI")

joindf <- joindf %>%
  mutate(
    bmsum = map(
      modx,
      reframe,
      across(
        any_of(contvars),
        function(x) {
          c(
            quantile(
              x,
              probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
              na.rm = TRUE
            ),
            sum(is.na(x))
          )
        }
      )
    )
  )

# Categorical variables, change accordingly
catvars <- covars[covars != "TSI"]

joindf <- joindf %>%
  mutate(
    covarsum = map(
      modx,
      ~.x |>
        count(across(any_of(catvars))) |>
        pivot_longer(
          -n,
          names_to = "covar_nm",
          values_to = "covar_val"
        ) |>
        group_by(covar_nm, covar_val) |>
        summarise(nc = sum(n), .groups = "drop")
    )
  )


resultdf <- joindf |>
  select(
    sex, reftime,
    LongDatSum,
    LongModSum,
    TotalN,
    bmsum,
    covarsum,
    modres
  )

# Cohort - Change accordingly
cohort <- "UKB"

# Saving submission file
saveRDS(
  resultdf,
  file = file.path(
    projfld, "data", cohort,
    paste0("resultdf_", cohort, ".rds")
  )
)

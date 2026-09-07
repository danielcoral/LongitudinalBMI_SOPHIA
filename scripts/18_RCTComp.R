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
library(ggplot2)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA/"

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
  "Menopause",
  # Ancestry
  "ancestry"
)

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
    )
  ) |>
  unnest(modlist)
print(regresdat)

sqrt(93.7 / 33.77)

# DPP BMI:
dpp_bmi <- 33.77
# DPP weight:
dpp_weight <- 93.7
# Average height in meters:
dpp_height <- sqrt(dpp_weight / dpp_bmi)
# Average weight loss at 2 years in kg:
dpp_wl_2y <- 5.39
wlpy <- dpp_wl_2y / 2
# Average weight loss in BMI units per year:
dpp_wloss_2yp <- round(
  (wlpy / (dpp_height^2)),
  2
)
dpp_wloss_2yp

dpptest <- regresdat |>
  filter(trait %in% c("GLU", "SBP", "TG", "ALT")) |>
  mutate(
    ATE = map2(
      mod1, bmdf,
      ~marginaleffects::avg_comparisons(
        .x,
        newdata = .y,
        variables = list(
          Slope = mutate(
            .y,
            high = Slope,
            low = pmin(-0.97, Slope),
            .keep = "none"
          )
        )
      )
    )
  )

dppres <- dpptest |>
  select(sex, trait, ATE) |>
  unnest(ATE) |>
  mutate(
    across(
      c(estimate, conf.low, conf.high),
      ~case_when(
        # Convert glucose from mmol/L to mg/dL
        trait == "GLU" ~ .x * 18,
        # Convert TG from mmol/L to mg/dL
        trait == "TG" ~ .x * 88.57,
        TRUE ~ .x
      )
    ),
    # Convert to per kg weight loss
    across(
      c(estimate, conf.low, conf.high),
      ~-.x / wlpy
    )
  )
dppres

dpp_eff <- tibble(
  term = "DPP",
  trait = c("GLU", "SBP", "TG", "ALT"),
  estimate = c(-0.57, -0.30, -2.27, -1),
  conf.low = c(-0.66, -0.40, -2.76, NA),
  conf.high = c(-0.48, -0.21, -1.78, NA)
)

options(repr.plot.width = 5, repr.plot.height = 2)
dppres |>
  select(term, sex, trait, estimate, conf.low, conf.high) |>
  mutate(term = paste("UKB", sex, sep = "_")) |>
  bind_rows(dpp_eff) |>
  ggplot(aes(estimate, term)) +
  geom_vline(xintercept = 0, linetype = "11") +
  geom_linerange(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      group = sex
    )
  ) +
  geom_point() +
  facet_grid(~trait, scales = "free_x") +
  labs(
    x = "Change in trait per kg weight loss (95% CI)",
    y = ""
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

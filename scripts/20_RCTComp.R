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

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA/"

agebmi <- c("age", "bmi")

bmrks <- c(
  # Biochemistry markers
  "ALT", "AST", "CRP", "DBP", "GGT", "HDL",
  "LDL", "GLU", "SBP", "SCR", "TG",
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

write_csv(
  tibble(baseline_bmi = dpp_bmi, baseline_weight_kg = dpp_weight,
         loss_two_years_kg = dpp_wl_2y,
         slope_bmi_per_year = -dpp_wloss_2yp),
  file.path(projfld, "data", "UKB", "results", "dpp_parameters.csv")
)

dpptest <- regresdat |>
  filter(trait %in% c("GLU", "SBP", "TG", "ALT")) |>
  mutate(
    ATE = map2(
      mod1, bmdf,
      function(model, data) {
        # Identity-link average contrast: compute its covariance directly.
        # Chunk design matrices to avoid large individual-gradient objects.
        stopifnot(model$family$link == "identity")
        contrast <- rep(0, length(coef(model)))
        chunks <- split(seq_len(nrow(data)),
                        ceiling(seq_len(nrow(data)) / 2000))
        for (rows in chunks) {
          high <- data[rows, ]
          low <- high
          low$Slope <- pmin(-dpp_wloss_2yp, low$Slope)
          contrast <- contrast + colSums(
            predict(model, high, type = "lpmatrix") -
              predict(model, low, type = "lpmatrix")
          ) / nrow(data)
        }
        estimate <- sum(contrast * coef(model))
        se <- sqrt(as.numeric(contrast %*% model$Vp %*% contrast))
        tibble(term = "Slope", estimate,
               conf.low = estimate - qnorm(0.975) * se,
               conf.high = estimate + qnorm(0.975) * se)
      }
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
    ),
    # Multiplication by a negative number reverses interval endpoints.
    lower = conf.high,
    conf.high = conf.low,
    conf.low = lower
  ) |>
  select(-lower)

dpp_eff <- tibble(
  term = "DPP",
  trait = c("GLU", "SBP", "TG", "ALT"),
  estimate = c(-0.57, -0.30, -2.27, -1),
  conf.low = c(-0.66, -0.40, -2.76, NA),
  conf.high = c(-0.48, -0.21, -1.78, NA)
)

dpp_plot <- dppres |>
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

saveRDS(
  dpp_plot,
  file.path(projfld, "data", "UKB", "processed_data", "plot_dpp_comparison.rds")
)
write_csv(
  dppres,
  file.path(projfld, "data", "UKB", "results", "table_dpp_comparison.csv")
)
write_lines(
  c(
    "# UKB projected biomarker change per kg weight loss",
    capture.output(print(dppres, row.names = FALSE))
  ),
  file.path(projfld, "data", "UKB", "results", "results_dpp_comparison.txt")
)

# Bootstrap

## Libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

## Working directory
projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

## Data
slopebmregdf <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf.rds")
)

## Covariates
covars <- c(
  # Smoking status
  "CurrSmoke",
  # Diabetes and insulin treatment
  "Diabetes", "Insulin",
  # Lipid-lowering treatment
  "LipidLower",
  # Antihypertensive treatment
  "HTMed",
  # Socioeconomic status
  "TSI", "CollegeEd",
  # Menopause
  "Menopause",
  # Ancestry
  "ancestry"

)

## Smooth terms
smterms <- c(
  "s(age, bs = \"cr\")",   # Non-linear effect of age
  "s(bmi, bs = \"cr\")",   # Non-linear effect of BMI
  "ti(age, bmi)",    # Variation in the effect of BMI due to age
  "s(Slope, bs = \"cr\")", # Non-linear effect of the BMI slope
  "ti(bmi, Slope)",  # Variation in the effect of BMI due to slope
  "ti(Slope, age)"   # Variation in the effect of the slope due to age
)

## Fitting models
slopebmregdf |>
  mutate(
    BTSTRP = pmap_chr(
      list(sex, trait, bmdf),
      function(sx, trt, x) {
        # Fitting models
        # Number of individuals
        ni <- nrow(x)
        # 1000 bootstraps
        modt <- tibble(B_ind = 1:1000) |>
          mutate(
            # Train set indices
            B_train = map(
              B_ind,
              function(i) {
                set.seed(i)
                sample(ni, round(ni * 0.8))
              }
            ),
            # Test set indices
            B_test = map(
              B_train,
              function(i) {
                (1:ni)[-i]
              }
            ),
            # Train and test datasets
            across(
              c(B_train, B_test),
              function(iveccol) {
                map(
                  iveccol,
                  function(ivec) {
                    x[ivec, ]
                  }
                )
              }
            ),
            # Fitting models
            B_mod = map(
              B_train,
              function(bd) {
                mgcv::bam(
                  reformulate(
                    termlabels = c(covars, smterms),
                    response = "value"
                  ),
                  data = bd,
                  drop.unused.levels = TRUE
                )
              }
            )
          )

        # Extracting mean squared error
        modt |>
          transmute(
            B_ind,
            B_mse1 = map_dbl(
              B_mod,
              function(mod) {
                mean(mod$residuals^2)
              }
            ),
            B_mse2 = map2_dbl(
              B_mod, B_test,
              function(mod, newd) {
                yhat <- predict(mod, newdata = newd)
                mean((newd$value - yhat)^2)
              }
            )
          ) |>
          write_tsv(
            file.path(
              projfld, "data", "UKB", "processed_data", "bootstrap_rec",
              paste0("MSEDAT_", sx, "_", trt, ".tsv")
            )
          )

        # Comparison between rapid gain vs rapid loss
        expectdf <- expand_grid(
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
        modt |>
          transmute(
            B_ind,
            modpred = map(
              B_mod,
              marginaleffects::comparisons,
              newdata = expectdf, vcov = FALSE,
              variables = list(Slope = c(-0.5, 0.5))
            )
          ) |>
          unnest(modpred) |>
          write_tsv(
            file.path(
              projfld, "data", "UKB", "processed_data", "bootstrap_rec",
              paste0("PREDS_", sx, "_", trt, ".tsv")
            )
          )

        # End!
        "Done"
      }
    )
  )
#End

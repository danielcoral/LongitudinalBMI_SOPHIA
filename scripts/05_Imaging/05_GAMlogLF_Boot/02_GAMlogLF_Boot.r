library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(mgcv)
library(patchwork)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

slopebmregdf_2ilf <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf_2ilf.rds")
)

print(slopebmregdf_2ilf)

bootdf <- slopebmregdf_2ilf |>
  transmute(
    sex, trait,
    BTSTRP = map2(
      sex, trait,
      function(sx, trt) {
        read_tsv(
          file.path(
            projfld, "data", "UKB", "processed_data", "bootstrap_2i",
            paste0("PREDS_", sx, "_", trt, ".tsv")
          ),
          show_col_types = FALSE
        )
      }
    )
  ) |>
  unnest(BTSTRP)

head(bootdf)

expectdiffdf <- expand_grid(
  CurrSmoke = 0.08,
  Diabetes = 0.05,
  LipidLower = 0.16,
  HTMed = 0.25,
  Insulin = 0.5,
  TSI = -1.3,
  CollegeEd = 0.3,
  Menopause = 0.4,
  age = c(45, 55, 65),
  bmi = c(22.5, 27.5, 32.5, 37.5),
  Slope = NA
)

pdiffdf <- slopebmregdf_2ilf |>
  transmute(
    sex, trait,
    modpred = map(
      mod1,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    thresh = map_dbl(bmdf, ~0.1 * sd(.x$value))
  ) |>
  unnest(modpred)

head(pdiffdf)

bmrks <- c(
  "LiverF_perc", "logLiverF_perc"
)
traitlabs <- gsub("_perc", "%", bmrks)

options(repr.plot.width = 3.5, repr.plot.height = 2)
bind_rows(
  Subsample = bootdf,
  Original = pdiffdf,
  .id = "Estimate"
) |>
  mutate(
    Estimate = factor(
      Estimate,
      levels = c("Original", "Subsample")
    ),
    trait = factor(trait, levels = bmrks, labels = traitlabs)
  ) |>
          ggplot(aes(bmi, estimate)) +
          geom_rect(
            aes(
              ymin = -thresh,
              ymax = thresh
            ),
            xmin = -Inf, xmax = Inf,
            fill = "grey70", alpha = 0.3,
            na.rm = TRUE
          ) +
          geom_hline(yintercept = 0, lty = "11", linewidth = .2) +
          geom_hline(
            aes(yintercept = -thresh),
            lty = "11", linewidth = .2,
            na.rm = TRUE
          ) +
          geom_hline(
            aes(yintercept = thresh),
            lty = "11", linewidth = .2,
            na.rm = TRUE
          ) +
          geom_point(
            aes(color = Estimate),
            size = 0.5
          ) +
          geom_linerange(
            aes(ymin = conf.low, ymax = conf.high),
            linewidth = .3, na.rm = TRUE
          ) +
          scale_color_manual(
            values = c("Original" = "black", "Subsample" = "#fec3ff")
          ) +
          ggh4x::facet_nested(
            trait ~ sex + age,
            scales = "free",
            strip = ggh4x::strip_nested(size = "variable")
          ) +
          labs(
            x = "BMI",
            y = "RG - RL difference",
            color = "Estimate"
          ) +
          theme_bw() +

  theme(
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2),
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = -10, b = -10),
    legend.key.width = unit(1, "mm")
  )

bootmse <- slopebmregdf_2ilf |>
  transmute(
    sex, trait,
    BTSTRP = map2(
      sex, trait,
      function(sx, trt){
        read_tsv(
          file.path(
            projfld, "data", "UKB", "processed_data", "bootstrap_2i",
            paste0("MSEDAT_", sx, "_", trt, ".tsv")
          ),
          show_col_types = FALSE
        )
      }
    )
  ) |>
  unnest(BTSTRP) |>
  mutate(
    across(c(B_mse1, B_mse2), sqrt),
    msediff = B_mse2 / B_mse1
  ) |>
  group_by(sex, trait) |>
  summarise(
    qmin = quantile(msediff, probs = .025),
    q50 = quantile(msediff, probs = .5),
    qmax = quantile(msediff, probs = .975),
    .groups = "drop"
  )
head(bootmse)

options(repr.plot.width = 3.5, repr.plot.height = 2)
bootmse |>
    ggplot(aes(q50, trait)) +
    geom_vline(xintercept = 1, lty = "11") +
    geom_vline(xintercept = c(.95, 1.05), lty = "11", color = "red") +
    #coord_cartesian(xlim = c(0.9, 1.1)) +
    geom_linerange(aes(xmin = qmin, xmax = qmax)) +
    geom_point() +
    facet_grid(~sex) +
    theme_bw() +
    labs(x = "Test/train MSE ratio", y = NULL)

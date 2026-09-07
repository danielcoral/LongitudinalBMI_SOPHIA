library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(mgcv)
library(patchwork)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

lmmtabrec <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
)

print(lmmtabrec)

bm1i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bm1i.tsv"),
  show_col_types = FALSE
)

head(bm1i)

cv1i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cv1i.tsv"),
  show_col_types = FALSE
)

head(cv1i)

ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  show_col_types = FALSE
) |>
  mutate(ancestry = factor(ancestry))
head(ancestrytab)

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

bmrks <- c(
  # Liver
  "ALT", "AST", "GGT",
  # Inflammatory
  "CRP",
  # Blood pressure
  "SBP", "DBP",
  # Lipids
  "HDL", "LDL", "TG",
  # Glycemic
  "GLU",
  # Renal
  "SCR",
  # Fat distribution
  "WHR"
)

slopebmregdf_delay <- bm1i |>
  # Joining biomarker and covariate data
  inner_join(cv1i, by = "eid") |>
  inner_join(ancestrytab, by = "eid") |>
  select(
    eid, sex, age, bmi,
    all_of(covars),
    all_of(bmrks)
  ) |>
  # Reshaping to long format
  pivot_longer(
    -c(eid, sex, age, bmi, all_of(covars)),
    names_to = "trait",
    values_drop_na = TRUE
  ) |>
  #  Nesting by sex
  nest(Dat = -sex) |>
  # Joining with trajectory models
  inner_join(lmmtabrec, by = "sex") |>
  # Joining with individual trajectory estimates
  transmute(
    sex,
    Dat = map2(indiv_est, Dat, inner_join, by = "eid"),
    Dat = map(Dat, nest, bmdf = -trait)
  ) |>
  unnest(Dat)

print(slopebmregdf_delay)

slopebmregdf_delay |>
  transmute(
    bmdf = map_dbl(bmdf, nrow)
  ) |>
  slice_max(bmdf)

smterms <- c(
  "s(age, bs = \"cr\")",   # Non-linear effect of age
  "s(bmi, bs = \"cr\")",   # Non-linear effect of BMI
  "ti(age, bmi)",        # Variation in the effect of BMI due to age
  "s(Slope, bs = \"cr\")", # Non-linear effect of the BMI slope
  "ti(bmi, Slope)",      # Variation in the effect of BMI due to slope
  "ti(Slope, age)"       # Variation in the effect of the slope due to age
)

slopebmregdf_delay <- slopebmregdf_delay |>
  mutate(
    mod1_delay = map(
      bmdf,
      ~bam(
        reformulate(
          termlabels = c(covars, smterms),
          response = "value"
        ),
        data = .x
      )
    )
  )

print(slopebmregdf_delay)

slopebmregdf <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf.rds")
)

expectdiffdf <- expand_grid(
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

grpdf <- bind_rows(
  tibble(grp1 = "Liver", trait = c("ALT", "AST", "GGT"), grp2 = 1),
  tibble(grp1 = "Lipids", trait = c("HDL", "LDL", "TG"), grp2 = 1),
  tibble(grp1 = "BP", trait = c("SBP", "DBP"), grp2 = 2),
  tibble(grp1 = "Glucose", trait = c("GLU"), grp2 = 2),
  tibble(grp1 = "Inflam", trait = c("CRP"), grp2 = 2),
  tibble(grp1 = "Kidney", trait = c("SCR"), grp2 = 2),
  tibble(grp1 = "Adipose", trait = c("WHR"), grp2 = 2)
)

pdiffdf <- slopebmregdf_delay |>
  select(sex, trait, mod1_delay) |>
  inner_join(slopebmregdf, by = c("sex", "trait")) |>
  transmute(
    sex, trait,
    modpred_all = map(
      mod1,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    modpred_delay = map(
      mod1_delay,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-0.5, 0.5))
    ),
    thresh = map_dbl(bmdf, ~0.1 * sd(.x$value))
  ) |>
  pivot_longer(
    -c(sex, trait, thresh),
    names_to = c(".value", "Subset"),
    names_sep = "_"
  )

options(repr.plot.width = 7, repr.plot.height = 4)
pdiffdf |>
  unnest(modpred) |>
  inner_join(grpdf, by = "trait") |>
  mutate(
    grp1 = factor(
      grp1,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney", "Adipose"
      )
    ),
    Association = factor(
      Subset,
      levels = c("all", "delay"),
      labels = c("Concurrent\n(Recruitment)", "Lagged\n(2 visit)")
    )
  ) |>
  nest(plotdf = -grp2) |>
  mutate(
    plotdf = map(
      plotdf,
      function(x) {
        x |>
          ggplot(aes(bmi, estimate)) +
          geom_rect(
            aes(
              ymin = -thresh,
              ymax = thresh
            ),
            xmin = -Inf, xmax = Inf,
            fill = "grey70", alpha = 0.3
          ) +
          geom_hline(yintercept = 0, lty = "11", linewidth = .2) +
          geom_hline(aes(yintercept = -thresh), lty = "11", linewidth = .2) +
          geom_hline(aes(yintercept = thresh), lty = "11", linewidth = .2) +
          geom_linerange(
            aes(ymin = conf.low, ymax = conf.high, group = Association),
            position = position_dodge(width = 1.25),
            linewidth = .3
          ) +
          geom_point(
            aes(group = Association, color = Association),
            position = position_dodge(width = 1.25),
            size = .5
          ) +
          ggh4x::facet_nested(
            grp1 + trait ~ sex + age,
            scales = "free",
            strip = ggh4x::strip_nested(size = "variable")
          ) +
          labs(
            x = "BMI",
            y = "RG - RL difference",
            color = "Association"
          ) +
          theme_bw()
      }
    ),
    plotdf = map2(
      grp2, plotdf,
      function(level, pl) {
        if (level == 2) {
          pl + theme(axis.title.y = element_blank())
        } else {
          pl
        }
      }
    )
  ) |>
  pull(plotdf) |>
  wrap_plots(nrow = 1, guides = "collect") &
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
    legend.key.width = unit(1, "mm"),
  )

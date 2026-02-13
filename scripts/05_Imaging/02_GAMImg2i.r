library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(mgcv)
library(patchwork)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

lmmtab2i <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "lmmtab2i.rds")
)

print(lmmtab2i)

img2i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "img2i.tsv"),
  show_col_types = FALSE
)

head(img2i)

cv2i <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cv2i.tsv"),
  show_col_types = FALSE
)

head(cv2i)

ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  col_types = "nf"
)
head(ancestrytab)

slopebmregdf_2i <- img2i |>
  inner_join(cv2i, by = "eid") |>
  inner_join(ancestrytab, by = "eid")

head(slopebmregdf_2i)

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
  "BF_perc",
  "ArmsF_perc",  "LegsF_perc", "TrunkF_perc",
  "AndF_perc", "GynF_perc",
  "LiverF_perc"
)

slopebmregdf_2i <- slopebmregdf_2i |>
  select(
    eid, sex, age, bmi,
    all_of(covars),
    all_of(bmrks)
  )

slopebmregdf_2i <- slopebmregdf_2i |>
  pivot_longer(
    -c(eid, sex, age, bmi, all_of(covars)),
    names_to = "trait",
    values_drop_na = TRUE
  )

slopebmregdf_2i <- slopebmregdf_2i |>
  nest(Dat = -sex) |>
  inner_join(lmmtab2i, by = "sex") |>
  transmute(
    sex,
    Dat = map2(indiv_est, Dat, inner_join, by = "eid"),
    Dat = map(Dat, nest, bmdf = -trait)
  ) |>
  unnest(Dat)

print(slopebmregdf_2i)

smterms <- c(
  "s(age, bs = \"cr\")",   # Non-linear effect of age
  "s(bmi, bs = \"cr\")",   # Non-linear effect of BMI
  "ti(age, bmi)",        # Variation in the effect of BMI due to age
  "s(Slope, bs = \"cr\")", # Non-linear effect of the BMI slope
  "ti(bmi, Slope)",      # Variation in the effect of BMI due to slope
  "ti(Slope, age)"       # Variation in the effect of the slope due to age
)

slopebmregdf_2i <- slopebmregdf_2i |>
  mutate(
    mod1 = map(
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

print(slopebmregdf_2i)

saveRDS(
  slopebmregdf_2i,
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf_2i.rds")
)

varexpdf_2i <- slopebmregdf_2i |>
  transmute(
    sex, trait,
    veval = map_dbl(mod1, ~summary(.x)[["r.sq"]])
  )

head(varexpdf_2i)

round(100 * summary(varexpdf_2i$veval), 2)

varexpdf_2i |>
  filter(grepl("_perc", trait)) |>
  mutate(veval = 100 * veval) |>
  pull(veval) |>
  summary() |>
  round(2)

varexpdf_2i |>
  filter(trait == "LiverF_perc") |>
  mutate(veval = 100 * veval) |>
  pull(veval) |>
  summary() |>
  round(2)

(traitlabs <- gsub("_perc", "%", bmrks))

options(repr.plot.width = 6, repr.plot.height = 3)
varexp_plot <- varexpdf_2i |>
  mutate(trait = factor(trait, levels = bmrks, labels = traitlabs)) |>
  ggplot(aes(veval, reorder(trait, veval))) +
  geom_col() +
  facet_grid(~sex) +
  theme_bw() +
  scale_x_continuous(labels = \(x) x * 100) +
  labs(x = "Variance Explained (%)", y = NULL)
varexp_plot

addvaluedf_2i <- slopebmregdf_2i |>
  mutate(

    # Model with covariates, age and BMI - removing slope
    mod_cvagebmi = map2(
      bmdf, mod1,
      function(x, m, covars, smterms) {
        # New smooth terms
        new_smterms <- c(
          "s(age, bs = \"cr\")",
          "s(bmi, bs = \"cr\")",
          "ti(age, bmi)"
        )
        # Extracting smooth penalties
        spterms <- c("s(age)", "s(bmi)", "ti(age,bmi)1", "ti(age,bmi)2")
        spvals <- m$sp[spterms]
        # Fitting reduced model
        bam(
          reformulate(
            termlabels = c(covars, new_smterms),
            response = "value"
          ),
          data = x,
          sp = spvals
        )
      },
      covars = covars,
      smterms = smterms
    ),

    # Model with covariates and age - removing BMI
    mod_cvage = map2(
      bmdf, mod1,
      function(x, m, covars, smterms) {
        # New smooth terms
        new_smterms <- "s(age, bs = \"cr\")"
        # Extracting smooth penalties
        spterms <- "s(age)"
        spvals <- m$sp[spterms]
        # Fitting reduced model
        bam(
          reformulate(
            termlabels = c(covars, new_smterms),
            response = "value"
          ),
          data = x,
          sp = spvals
        )
      },
      covars = covars,
      smterms = smterms
    ),

    # Model with covariates alone
    mod_cv = map(
      bmdf,
      function(x, covars) {
        # Fitting reduced model
        bam(
          reformulate(
            termlabels = covars,
            response = "value"
          ),
          data = x
        )
      },
      covars = covars
    ),

    # Model with covariates, age BMI and slope, removing BMIxSlope
    mod_cvagebmislope = map2(
      bmdf, mod1,
      function(x, m, covars, smterms) {
        # New smooth terms
        new_smterms <- c(
          "s(age, bs = \"cr\")",
          "s(bmi, bs = \"cr\")",
          "ti(age, bmi)",
          "s(Slope, bs = \"cr\")",
          "ti(Slope, age)"
        )
        # Extracting smooth penalties
        spterms <- c(
          "s(age)", "s(bmi)",
          "ti(age,bmi)1", "ti(age,bmi)2",
          "s(Slope)",
          "ti(Slope,age)1", "ti(Slope,age)2"
        )
        spvals <- m$sp[spterms]
        # Fitting reduced model
        bam(
          reformulate(
            termlabels = c(covars, new_smterms),
            response = "value"
          ),
          data = x,
          sp = spvals
        )
      },
      covars = covars,
      smterms = smterms
    )

  )

print(addvaluedf_2i)

addedvalslope_2i <- addvaluedf_2i |>
  transmute(
    sex,
    trait,
    P_addslope = map2_dbl(
      mod_cvagebmi, mod1,
      ~anova(.x, .y)[2, 6]
    ),
    FDR = p.adjust(P_addslope, "fdr") < 0.05,
    FAI = map2_dbl(
      mod_cvagebmi, mod1,
      function(mr, mc) {
        llnull <- mc$null.deviance
        llc <- deviance(mc) - llnull
        llr <- deviance(mr) - llnull
        round(100 * (1 - (llr / llc)), 2)
      }
    )
  )

addedvalslope_2i

write_tsv(
  addedvalslope_2i,
  file.path(projfld, "data", "UKB", "processed_data", "addedvalslope_2i.tsv")
)

addedvalslopeint_2i <- addvaluedf_2i |>
  transmute(
    sex,
    trait,
    P_addslopeint = map2_dbl(
      mod_cvagebmislope, mod1,
      ~anova(.x, .y)[2, 6]
    ),
    FDR = p.adjust(P_addslopeint, "fdr") < 0.05,
    FAI = map2_dbl(
      mod_cvagebmislope, mod1,
      function(mr, mc) {
        llnull <- mc$null.deviance
        llc <- deviance(mc) - llnull
        llr <- deviance(mr) - llnull
        round(100 * (1 - (llr / llc)), 2)
      }
    )
  )

addedvalslopeint_2i

FAIdf_2i <- addvaluedf_2i |>
  transmute(
    sex,
    trait,
    across(
      c(mod_cv, mod_cvage, mod_cvagebmi),
      function(m2) {
        map2_dbl(
          m2, mod1,
          function(mr, mc) {
            llnull <- mc$null.deviance
            llc <- deviance(mc) - llnull
            llr <- deviance(mr) - llnull
            round(100 * llr / llc, 2)
          }
        )
      }
    ),
    FAI_cv = mod_cv,
    FAI_age = mod_cvage - mod_cv,
    FAI_bmi = mod_cvagebmi - mod_cvage,
    FAI_slope = 100 - mod_cvagebmi
  ) |>
  select(-starts_with("mod"))

head(FAIdf_2i)

summary(FAIdf_2i$FAI_bmi)

FAIdf_2i |>
  filter(grepl("_perc", trait)) |>
  pull(FAI_bmi) |>
  summary()

FAIdf_2i |>
  filter(trait == "LiverF_perc") |>
  pull(FAI_bmi) |>
  summary()

slice_max(FAIdf_2i, FAI_slope, by = "sex")

options(repr.plot.width = 6, repr.plot.height = 3)
fai_plot <- FAIdf_2i |>
  pivot_longer(
    starts_with("FAI_"),
    names_to = c(".value", "term"),
    names_sep = "_"
  ) |>
  mutate(
    term = case_match(
      term,
      "cv" ~ "Covariates",
      "age" ~ "Age",
      "bmi" ~ "BMI",
      "slope" ~ "Slope"
    ),
    term = factor(
      term,
      levels = rev(c("Covariates", "Age", "BMI", "Slope"))
    ),
    trait = gsub("_perc", "%", trait)
  ) |>
  group_by(trait) |>
  mutate(slopemax = max(FAI[term == "Slope"])) |>
  ungroup() |>
  ggplot(aes(reorder(trait, slopemax), FAI)) +
  geom_col(aes(fill = term)) +
  facet_grid(~sex) +
  theme_bw() +
  labs(y = "Proportion of explained variance", x = NULL) +
  guides(fill = guide_legend(title = NULL, reverse = TRUE)) +
  coord_flip()
fai_plot

bmivslope_2i <- FAIdf_2i |>
  transmute(
    sex, trait,
    prop_bmi = round(100 * FAI_bmi / (FAI_bmi + FAI_slope), 2),
    prop_slope = round(100 * FAI_slope / (FAI_bmi + FAI_slope), 2)
  )


bmivslope_2i

options(repr.plot.width = 6, repr.plot.height = 3)
relfai_plot <- bmivslope_2i |>
  pivot_longer(
    starts_with("prop_"),
    names_to = c(".value", "term"),
    names_sep = "_"
  ) |>
  mutate(
    term = case_match(
      term,
      "bmi" ~ "BMI",
      "slope" ~ "Slope"
    ),
    term = factor(
      term,
      levels = rev(c("BMI", "Slope"))
    ),
    trait = gsub("_perc", "%", trait)
  ) |>
  group_by(trait) |>
  mutate(propmax = max(prop[term == "Slope"])) |>
  ungroup() |>
  ggplot(aes(reorder(trait, propmax), prop)) +
  geom_col(aes(fill = term)) +
  scale_fill_manual(values = scales::hue_pal()(4)[1:2]) +
  facet_grid(~sex) +
  theme_bw() +
  labs(y = "Relative proportion of explained variance", x = NULL) +
  guides(fill = guide_legend(title = NULL, reverse = TRUE)) +
  coord_flip()
relfai_plot

expectdf <- expand_grid(
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
  Slope = c(-0.5, 0, 0.5),
  ancestry = "EUR"
)

predtab <- slopebmregdf_2i |>
  transmute(
    sex, trait,
    mod1_pred = map(
      mod1,
      ~gratia::fitted_values(.x, data = expectdf)
    )
  ) |>
  unnest(mod1_pred)

head(predtab)

options(repr.plot.width = 7, repr.plot.height = 4)
predtab |>
  mutate(trait = factor(trait, levels = bmrks)) |>
  ggplot(aes(bmi, .fitted)) +
  geom_linerange(
    aes(ymin = .lower_ci, ymax = .upper_ci, group = Slope),
    position = position_dodge(width = 1.25),
    linewidth = .3
  ) +
  geom_point(
    aes(group = Slope, color = factor(Slope)),
    position = position_dodge(width = 1.25),
    size = .5
  ) +
  ggh4x::facet_nested(
    trait ~ sex + age,
    scales = "free",
    strip = ggh4x::strip_nested(size = "variable")
  ) +
  labs(
    x = "BMI",
    y = "Expected trait value",
    color = "Slope (kg/m2/yr)"
  ) +
  theme_bw() +
  theme(
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = -10, b = -10),
    legend.key.width = unit(1, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

expectdiffdf <- expectdf |>
  mutate(Slope = NA) |>
  unique()

pdiffdf_2i <- slopebmregdf_2i |>
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

head(pdiffdf_2i)

options(repr.plot.width = 3.5, repr.plot.height = 3)
pdiff_plot <- pdiffdf_2i |>
  mutate(trait = factor(trait, levels = bmrks, labels = traitlabs)) |>
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
    aes(ymin = conf.low, ymax = conf.high),
    linewidth = .3
  ) +
  geom_point(size = .5) +
  ggh4x::facet_nested(
    trait ~ sex + age,
    scales = "free",
    strip = ggh4x::strip_nested(size = "variable")
  ) +
  labs(
    x = "BMI",
    y = "RG - RL difference"
  ) +
  theme_bw() +
  theme(
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )
pdiff_plot

options(repr.plot.width = 7, repr.plot.height = 5)
design <- "
ADDD
BDDD
CDDD
"
patchwork::wrap_plots(
  A = varexp_plot + theme(axis.text = element_text(size = 5)),
  B = fai_plot + theme(
    legend.key.size = unit(2, "mm"),
    legend.position = "top",
    legend.justification = "left"
  ) + guides(
    fill = guide_legend(title = "Model\ncomponents", nrow = 2, reverse = TRUE)
  ),
  C = relfai_plot + theme(legend.position = "none"),
  D = pdiff_plot,
  design = design, guides = "keep"
) &
  patchwork::plot_annotation(tag_levels = list("A", "")) &
  theme(
    plot.tag = element_text(size = 5),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 4),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    legend.margin = margin(t = 0, b = 0),
    strip.text = element_text(size = 5, margin = margin(1, 1, 1, 1))
  )

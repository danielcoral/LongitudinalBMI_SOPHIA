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

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)

head(bmrec)

cvrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv"),
  show_col_types = FALSE
)

head(cvrec)

ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  show_col_types = FALSE
) |>
  mutate(ancestry = factor(ancestry))
head(ancestrytab)

slopebmregdf <- bmrec |>
  inner_join(cvrec, by = "eid") |>
  inner_join(ancestrytab, by = "eid")

head(slopebmregdf)

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

slopebmregdf <- slopebmregdf |>
  select(
    eid, sex, age, bmi,
    all_of(covars),
    all_of(bmrks)
  )

slopebmregdf <- slopebmregdf |>
  pivot_longer(
    -c(eid, sex, age, bmi, all_of(covars)),
    names_to = "trait",
    values_drop_na = TRUE
  )

slopebmregdf <- slopebmregdf |>
  nest(Dat = -sex) |>
  inner_join(lmmtabrec, by = "sex") |>
  transmute(
    sex,
    Dat = map2(indiv_est, Dat, inner_join, by = "eid"),
    Dat = map(Dat, nest, bmdf = -trait)
  ) |>
  unnest(Dat)

print(slopebmregdf)

smterms <- c(
  "s(age, bs = \"cr\")",   # Non-linear effect of age
  "s(bmi, bs = \"cr\")",   # Non-linear effect of BMI
  "ti(age, bmi)",        # Variation in the effect of BMI due to age
  "s(Slope, bs = \"cr\")", # Non-linear effect of the BMI slope
  "ti(bmi, Slope)",      # Variation in the effect of BMI due to slope
  "ti(Slope, age)"       # Variation in the effect of the slope due to age
)

slopebmregdf <- slopebmregdf |>
  mutate(
    mod1 = map(
      bmdf,
      ~bam(
        reformulate(
          termlabels = c(covars, smterms),
          response = "value"
        ),
        data = .x,
        drop.unused.levels = TRUE
      )
    )
  )

print(slopebmregdf)

saveRDS(
  slopebmregdf,
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf.rds")
)

slopebmregdf <- readRDS(
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregdf.rds")
)

slopebmregdf |>
  transmute(
    sex, trait,
    bmdf = map(bmdf, count, ancestry),
    bmdf = map(bmdf, mutate, prop = round(100 * n / sum(n), 1))
  ) |>
  unnest(bmdf) |>
  group_by(ancestry) |>
  slice_max(prop, n = 1, with_ties = FALSE)

slopebmregres <- slopebmregdf |>
  transmute(
    sex, trait,
    mod1 = map(mod1, broom::tidy)
  ) |>
  unnest(mod1)

head(slopebmregres)

write_tsv(
  slopebmregres,
  file.path(projfld, "data", "UKB", "processed_data", "slopebmregres.tsv")
)

modelkcheck <- slopebmregdf |>
  transmute(
    sex, trait,
    Check_K = map(mod1, mgcv::k.check),
    Check_K = map(Check_K, data.frame),
    Check_K = map(Check_K, tibble::rownames_to_column, var = "term")
  ) |>
  unnest(Check_K) |>
  select(-p.value)

head(modelkcheck)

write_tsv(
  modelkcheck,
  file.path(projfld, "data", "UKB", "processed_data", "modelkcheck.tsv")
)

summary(modelkcheck$k.index)

varexpdf <- slopebmregdf |>
  transmute(
    sex, trait,
    veval = map_dbl(mod1, ~summary(.x)[["r.sq"]])
  )

head(varexpdf)

write_tsv(
  varexpdf,
  file.path(
    projfld, "data", "UKB", "processed_data", "varexpdf.tsv"
  )
)

round(100 * summary(varexpdf$veval), 2)

options(repr.plot.width = 6, repr.plot.height = 3)
varex_plot <- varexpdf |>
  ggplot(aes(veval, reorder(trait, veval))) +
  geom_col() +
  facet_grid(~sex) +
  theme_bw() +
  scale_x_continuous(labels = \(x) x * 100) +
  labs(x = "Variance Explained (%)", y = NULL)
varex_plot

addvaluedf <- slopebmregdf |>
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

print(addvaluedf)

addedvalslope <- addvaluedf |>
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

head(addedvalslope)

write_tsv(
  addedvalslope,
  file.path(projfld, "data", "UKB", "processed_data", "addedvalslope.tsv")
)

addedvalslopeint <- addvaluedf |>
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

addedvalslopeint

write_tsv(
  addedvalslopeint,
  file.path(projfld, "data", "UKB", "processed_data", "addedvalslopeint.tsv")
)

FAIdf <- addvaluedf |>
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

head(FAIdf)

options(repr.plot.width = 6, repr.plot.height = 3)
fai_plot <- FAIdf |>
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
    )
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

bmivslope <- FAIdf |>
  transmute(
    sex, trait,
    prop_bmi = round(100 * FAI_bmi / (FAI_bmi + FAI_slope), 2),
    prop_slope = round(100 * FAI_slope / (FAI_bmi + FAI_slope), 2)
  )


bmivslope

options(repr.plot.width = 6, repr.plot.height = 3)
relfai_plot <- bmivslope |>
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
    )
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

predtab <- slopebmregdf |>
  transmute(
    sex, trait,
    mod1_pred = map(
      mod1,
      ~gratia::fitted_values(.x, data = expectdf)
    )
  ) |>
  unnest(mod1_pred)

head(predtab)

grpdf <- bind_rows(
  tibble(grp1 = "Liver", trait = c("ALT", "AST", "GGT"), grp2 = 1),
  tibble(grp1 = "Lipids", trait = c("HDL", "LDL", "TG"), grp2 = 1),
  tibble(grp1 = "BP", trait = c("SBP", "DBP"), grp2 = 2),
  tibble(grp1 = "Glucose", trait = c("GLU"), grp2 = 2),
  tibble(grp1 = "Inflam", trait = c("CRP"), grp2 = 2),
  tibble(grp1 = "Kidney", trait = c("SCR"), grp2 = 2),
  tibble(grp1 = "Adipose", trait = c("WHR"), grp2 = 2)
)

options(repr.plot.width = 7, repr.plot.height = 4)
predtab |>
  inner_join(grpdf, by = "trait") |>
  mutate(
    grp1 = factor(
      grp1,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney", "Adipose"
      )
    )
  ) |>
  nest(plotdf = -grp2) |>
  mutate(
    plotdf = map(
      plotdf,
      function(x) {
        x |>
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
            grp1 + trait ~ sex + age,
            scales = "free",
            strip = ggh4x::strip_nested(size = "variable")
          ) +
          labs(
            x = "BMI",
            y = "Expected trait value",
            color = "Slope (kg/m2/yr)"
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

pdiffdf <- slopebmregdf |>
  transmute(
    sex, trait,
    modpred = map(
      mod1,
      marginaleffects::comparisons,
      newdata = expectdiffdf,
      variables = list(Slope = c(-.5, .5))
    ),
    thresh = map_dbl(bmdf, ~0.1 * sd(.x$value))
  ) |>
  unnest(modpred)

head(pdiffdf)

options(repr.plot.width = 7, repr.plot.height = 4)
pdiff_plot <- pdiffdf |>
  inner_join(grpdf, by = "trait") |>
  mutate(
    grp1 = factor(
      grp1,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose", "Inflam", "Kidney", "Adipose"
      )
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
            aes(ymin = conf.low, ymax = conf.high),
            linewidth = .3
          ) +
          geom_point(size = .5) +
          ggh4x::facet_nested(
            grp1 + trait ~ sex + age,
            scales = "free",
            strip = ggh4x::strip_nested(size = "variable")
          ) +
          labs(
            x = "BMI",
            y = "RG - RL difference"
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
  wrap_plots(nrow = 1) &
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

options(repr.plot.width = 3, repr.plot.height = 6)
bivardistrib_plot <- lmmtabrec |>
  transmute(
    sex,
    cordat = map(LongD, select, eid, bmi0),
    cordat = map(cordat, unique),
    cordat = map2(cordat, indiv_est, inner_join, by = "eid")
  ) |>
  unnest(cordat) |>
  ggplot(aes(x = bmi0, y = Slope)) +
  ggdensity::geom_hdr() +
  geom_smooth(method = "lm", formula = "y ~ x", color = "red", linewidth = .5) +
  geom_vline(
    xintercept = c(22.5, 27.5, 32.5, 37.5, 42.5),
    linetype = "11", color = "magenta",
    linewidth = .5
  ) +
  geom_hline(
    yintercept = c(-.5, 0, .5), linetype = "11",
    color = "darkgreen", linewidth = .5
  ) +
  geom_point(
    data = expand_grid(
      bmi0 = c(22.5, 27.5, 32.5, 37.5),
      Slope = c(-.5, 0, .5)
    ),
    color = "blue", size = .5
  ) +
  facet_wrap(~sex, ncol = 1) +
  theme_bw() +
  labs(
    x = "BMI at recruitment (kg/m²)",
    y = "Individual slope (kg/m² per year)",
    alpha = "Probability\ndensity"
  ) +
  theme(
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "top",
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  coord_cartesian(xlim = c(20, 45), ylim = c(-2.5, 2.5))
bivardistrib_plot

options(repr.plot.width = 7, repr.plot.height = 5)
design <- "
ACCDD
AEEEE
BEEEE
BEEEE
"
patchwork::wrap_plots(
  A = bivardistrib_plot + theme(
    legend.position = "top",
    legend.key.size = unit(2, "mm"),
    legend.text.position = "top"
  ) +
    guides(
      alpha = guide_legend(reverse = TRUE)
    ),
  B = varex_plot + theme(axis.text = element_text(size = 5)),
  C = fai_plot + theme(
    legend.key.size = unit(2, "mm"),
    legend.position = "top",
    legend.justification = "left"
  ) + guides(
    fill = guide_legend(title = "Model\ncomponents", nrow = 2, reverse = TRUE)
  ),
  D = relfai_plot + theme(legend.position = "none"),
  E = pdiff_plot,
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

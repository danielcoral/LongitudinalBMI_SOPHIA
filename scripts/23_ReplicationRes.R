# Cross-cohort descriptive summaries from the saved cohort result files.
# Inputs: data/{ALSPAC,RS,UKB}/resultdf_*.rds.
# Outputs: plots/cross_cohort_biomarker_distributions.png and
#   plots/cross_cohort_comorbidities.png; table_cross_cohort_sample_sizes.csv
#   and table_cross_cohort_biomarker_distributions.csv in data/UKB/results/.
# Slope contrasts across cohorts are in 37_WithinCohortContrasts.R.
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

restab <- expand_grid(
  cohort = c("ALSPAC", "RS", "UKB")
) |>
  mutate(
    resdf = map(
      cohort,
      ~readRDS(
        file.path(
          projfld, "data", .x,
          paste0("resultdf_", .x, ".rds")
        )
      )
    )
  )
print(restab)

restabl <- restab |>
  unnest(resdf) |>
  mutate(modres = map(modres, ~filter(.x, trait != "WHR")))

print(restabl)

totalss <- restabl |>
  select(cohort, reftime, TotalN) |>
  group_by(cohort, reftime) |>
  summarise(n = sum(TotalN), .groups = "drop")
totalss

bms <- restabl |>
  transmute(
    cohort = ifelse(
      cohort == "UKB" & reftime == "2i",
      paste(cohort, reftime, sep = "_"),
      cohort
    ),
    sex, bmsum
  ) |>
  mutate(
    bmsum = map(
      bmsum, mutate,
      valuenm = c(
        "qmin", "q25", "q50", "q75", "qmax", "nmiss"
      )
    )
  ) |>
  unnest(bmsum) |>
  pivot_longer(
    -c(cohort, sex, valuenm),
    names_to = "trait",
    values_to = "traitval",
    values_drop_na = TRUE,
    cols_vary = "slowest"
  )
print(bms)

bms |>
  filter(trait %in% c("age", "bmi")) |>
  pivot_wider(
    names_from = c(trait, valuenm),
    values_from = traitval
  )

bms |>
  filter(trait == "Slope") |>
  select(-trait) |>
  pivot_wider(
    names_from = valuenm,
    values_from = traitval
  ) |>
  mutate(across(where(is.numeric), ~round(.x, 2)))

bmrks <- c(
  # Biochemistry markers
  "ALT", "AST", "CRP", "DBP", "GGT", "HDL",
  "LDL", "GLU", "SBP", "SCR", "TG",
  # DEXA markers
  "BF_perc", "AndF_perc", "GynF_perc",
  "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
)

outbms <- bms |>
  filter(
    !(
      cohort == "ALSPAC" &
        grepl("F_perc$", trait) &
        !grepl("TF_perc$", trait)
    )
  ) |>
  mutate(
    cohort = gsub("_2i$", "", cohort),
    trait = gsub("TF_perc$", "F_perc", trait),
    trait = case_match(
      trait,
      "FG" ~ "GLU",
      "PCr_NMR" ~ "SCR",
      .default = trait
    )
  ) |>
  filter(trait %in% bmrks) |>
  pivot_wider(names_from = valuenm, values_from = traitval) |>
  # ALSPAC android fat percentage was provided inverted (as lean percentage).
  # Inverting to fat% (fat = 100 - lean) reverses the quantile order, so the
  # minimum and maximum swap and the 25th/75th percentiles swap.
  mutate(
    rev_min = 100 - qmax, rev_q25 = 100 - q75, rev_q50 = 100 - q50,
    rev_q75 = 100 - q25, rev_max = 100 - qmin
  ) |>
  mutate(
    qmin = if_else(cohort == "ALSPAC" & trait == "AndF_perc", rev_min, qmin),
    q25 = if_else(cohort == "ALSPAC" & trait == "AndF_perc", rev_q25, q25),
    q50 = if_else(cohort == "ALSPAC" & trait == "AndF_perc", rev_q50, q50),
    q75 = if_else(cohort == "ALSPAC" & trait == "AndF_perc", rev_q75, q75),
    qmax = if_else(cohort == "ALSPAC" & trait == "AndF_perc", rev_max, qmax)
  ) |>
  select(-rev_min, -rev_q25, -rev_q50, -rev_q75, -rev_max)

head(outbms)

grpdf <- bind_rows(
  tibble(grp1 = "Liver", trait = c("ALT", "AST", "GGT"), grp2 = 1),
  tibble(grp1 = "Lipids", trait = c("HDL", "LDL", "TG"), grp2 = 1),
  tibble(grp1 = "BP", trait = c("SBP", "DBP"), grp2 = 1),
  tibble(grp1 = "Glucose", trait = c("GLU"), grp2 = 1),
  tibble(grp1 = "Inflam", trait = c("CRP"), grp2 = 2),
  tibble(grp1 = "Kidney", trait = c("SCR"), grp2 = 2),
  tibble(
    grp1 = "Imaging",
    trait = c(
      "BF_perc", "AndF_perc", "GynF_perc",
      "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
    ),
    grp2 = 2
  )
)

options(repr.plot.width = 5, repr.plot.height = 4)
cohort_biomarker_distributions_plot <- outbms |>
  inner_join(grpdf, by = "trait") |>
  mutate(
    trait = gsub("_perc$", "%", trait),
    grp1 = factor(
      grp1,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose",
        "Inflam", "Kidney", "Adipose", "Imaging"
      )
    )
  ) |>
  nest(Dat = -grp2) |>
  mutate(
    Plt = map(
      Dat,
      ~.x |>
        ggplot(
          aes(cohort, q50)
        ) +
        geom_boxplot(
          aes(
            min = qmin, lower = q25,
            middle = q50,
            upper = q75, max = qmax
          ),
          stat = "identity"
        ) +
        ggh4x::facet_nested(
          grp1 + trait ~ sex, scales = "free_y",
          strip = ggh4x::strip_nested(size = "variable")
        ) +
        labs(x = "Cohort", y = "Biomarker value") +
        theme_bw()
    ) |>
      map2(
        grp2,
        ~ if (.y == 2) {
          .x +

            theme(
              axis.title.y = element_blank()
            )
        } else {
          .x
        }
      )
  ) |>
  pull(Plt) |>
  patchwork::wrap_plots(nrow = 1) &
  theme(
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

restabl |>
  filter(is.na(reftime) | reftime == "rec") |>
  transmute(
    cohort, sex,
    covarsum = map(
      covarsum, mutate,
      covar_val = as.character(covar_val)
    )
  ) |>
  unnest(covarsum) |>
  group_by(cohort, sex, covar_nm) |>
  mutate(prop = nc / sum(nc)) |>
  ungroup() |>
  ggplot(
    aes(cohort, prop)
  ) +
  geom_col(aes(group = covar_val, fill = covar_val), position = "stack") +
  facet_grid(sex ~ covar_nm, scales = "free") +
  labs(x = NULL, y = "Proportion", fill = "Value") +
  theme_bw() +
  theme(
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 4),
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

cohort_comorbidities_plot <- restabl |>
  filter(is.na(reftime) | reftime == "rec") |>
  transmute(
    cohort, sex,
    covarsum = map(
      covarsum, mutate,
      covar_val = as.character(covar_val)
    )
  ) |>
  unnest(covarsum) |>
  filter(
    covar_nm %in% c(
      "Menopause", "Insulin",
      "HTMed", "Diabetes", "CurrSmoke"
    ),
    !(covar_nm == "Menopause" & sex == "Male")
  ) |>
  group_by(cohort, sex, covar_nm) |>
  mutate(prop = nc / sum(nc)) |>
  ungroup() |>
  ggplot(
    aes(cohort, prop)
  ) +
  geom_col(aes(group = covar_val, fill = covar_val), position = "stack") +
  facet_grid(sex ~ covar_nm, scales = "free") +
  labs(x = NULL, y = "Proportion", fill = "Value") +
  theme_bw() +
  theme(
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 4),
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

cohort_biomarker_distributions_plot
ggsave(
  file.path(projfld, "plots", "cross_cohort_biomarker_distributions.png"),
  plot = cohort_biomarker_distributions_plot, width = 5, height = 4,
  units = "in", dpi = 300, bg = "white"
)
cohort_comorbidities_plot
ggsave(
  file.path(projfld, "plots", "cross_cohort_comorbidities.png"),
  plot = cohort_comorbidities_plot, width = 5, height = 4,
  units = "in", dpi = 300, bg = "white"
)
write_csv(
  totalss,
  file.path(
    projfld, "data", "UKB", "results", "table_cross_cohort_sample_sizes.csv"
  )
)
write_csv(
  outbms,
  file.path(
    projfld, "data", "UKB", "results",
    "table_cross_cohort_biomarker_distributions.csv"
  )
)

library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(mgcv)
library(metagam)

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
  unnest(resdf)

print(restabl)

totalss <- restabl |>
  select(cohort, reftime, TotalN) |>
  group_by(cohort, reftime) |>
  summarise(n = sum(TotalN), .groups = "drop")
totalss

lds <- restabl |>
  select(cohort, sex, reftime, LongDatSum) |>
  unnest(LongDatSum)
lds

lds |>
  transmute(
    cohort = ifelse(
      cohort == "UKB" & reftime == "2i",
      paste(cohort, reftime, sep = "_"),
      cohort
    ),
    sex,
    N = TotalInd
  ) |>
  group_by(cohort) |>
  mutate(prop = N / sum(N)) |>
  ungroup() |>
  arrange(cohort, sex)

lms <- restabl |>
  transmute(
    cohort = ifelse(
      cohort == "UKB" & reftime == "2i",
      paste(cohort, reftime, sep = "_"),
      cohort
    ),
    sex, LongModSum
  ) |>
  unnest(LongModSum) |>
  select(-std.error, -statistic) |>
  arrange(cohort, sex)
lms

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
  "LDL", "GLU", "SBP", "SCR", "TG", "WHR",
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
  pivot_wider(names_from = valuenm, values_from = traitval)

head(outbms)

grpdf <- bind_rows(
  tibble(grp1 = "Liver", trait = c("ALT", "AST", "GGT"), grp2 = 1),
  tibble(grp1 = "Lipids", trait = c("HDL", "LDL", "TG"), grp2 = 1),
  tibble(grp1 = "BP", trait = c("SBP", "DBP"), grp2 = 1),
  tibble(grp1 = "Glucose", trait = c("GLU"), grp2 = 1),
  tibble(grp1 = "Inflam", trait = c("CRP"), grp2 = 2),
  tibble(grp1 = "Kidney", trait = c("SCR"), grp2 = 2),
  tibble(grp1 = "Adipose", trait = c("WHR"), grp2 = 2),
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
outbms |>
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


predtab <- restabl |>
  mutate(
    cohort, sex, modres,
    # Median age
    agem = round(map_dbl(bmsum, ~.x[3, "age", drop = TRUE])),
    # Median and IQR of BMI
    bmi_med = round(map_dbl(bmsum, ~.x[3, "bmi", drop = TRUE]), 1),
    bmi_lwr = round(map_dbl(bmsum, ~.x[2, "bmi", drop = TRUE]), 1),
    bmi_upr = round(map_dbl(bmsum, ~.x[4, "bmi", drop = TRUE]), 1),
    # IQR of Slope
    Slope_lwr = round(map_dbl(bmsum, ~.x[2, "Slope", drop = TRUE]), 2),
    Slope_upr = round(map_dbl(bmsum, ~.x[4, "Slope", drop = TRUE]), 2),
    .keep = "none"
  ) |>
  # Choosing parameters for each cohort
  group_by(cohort) |>
  mutate(
    agem = min(agem),
    bmi_med = median(bmi_med),
    across(c(Slope_lwr, bmi_lwr), min),
    across(c(Slope_upr, bmi_upr), max)
  ) |>
  ungroup()
print(predtab)

predtab <- predtab |>
  unnest(modres)
print(predtab)

predtab <- predtab |>
  mutate(
    # Predictors for each model
    NewDat = map(mod1, pluck, "var.summary") |>
      map(data.frame) |>
      # Median values for continuous data
      map(slice, 2) |>
      # Adding median age
      map2(agem, ~mutate(.x, age = .y)),
    # Adding BMI IQRs
    NewDat = map(NewDat, select, -bmi),
    NewDat = pmap(
      list(NewDat, bmi_lwr, bmi_med, bmi_upr),
      ~expand_grid(..1, bmi = c(..2, ..3, ..4))
    )
  ) |>
  select(-c(agem, bmi_med, bmi_lwr, bmi_upr))
print(predtab)

predtab <- predtab |>
  mutate(
    SlopeExtDiffs = pmap(
      list(mod1, NewDat, Slope_upr, Slope_lwr),
      function(md, nd, sup, slw) {
        marginaleffects::comparisons(
          md,
          newdata = nd,
          variables = list(Slope = c(sup, slw))
        )
      }
    )
  )

predtab <- predtab |>
  mutate(
    cohort, sex, trait,
    Slope_upr, Slope_lwr,
    SlopeExtDiffs = map(
      SlopeExtDiffs,
      ~.x |>
        mutate(
          bmi, estimate, conf.low, conf.high,
          bmiq = c("25%", "50%", "75%"),
          .keep = "none"
        )
    ),
    .keep = "none"
  ) |>
  unnest(SlopeExtDiffs)
head(predtab)

outbmsdifs <- predtab |>
  # Harmonizing trait names before filtering
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
  # Selecting biomarkers of interest
  filter(trait %in% bmrks)

grpdf <- bind_rows(
  tibble(grp1 = "Liver", trait = c("ALT", "AST", "GGT"), grp2 = 1),
  tibble(grp1 = "Lipids", trait = c("HDL", "LDL", "TG"), grp2 = 1),
  tibble(grp1 = "BP", trait = c("SBP", "DBP"), grp2 = 1),
  tibble(grp1 = "Glucose", trait = c("GLU"), grp2 = 1),
  tibble(grp1 = "Inflam", trait = c("CRP"), grp2 = 2),
  tibble(grp1 = "Kidney", trait = c("SCR"), grp2 = 2),
  tibble(grp1 = "Adipose", trait = c("WHR"), grp2 = 2),
  tibble(
    grp1 = "Imaging",
    trait = c(
      "BF_perc", "AndF_perc", "GynF_perc",
      "ArmsF_perc", "LegsF_perc", "TrunkF_perc"
    ),
    grp2 = 2
  )
)

options(repr.plot.width = 5, repr.plot.height = 3.5)
outbmsdifs |>
  inner_join(grpdf, by = "trait") |>
  mutate(
    trait = gsub("_perc$", "%", trait),
    grp1 = factor(
      grp1,
      levels = c(
        "Liver", "Lipids", "BP", "Glucose",
        "Inflam", "Kidney", "Adipose", "Imaging"
      )
    ),
    cohort = factor(
      cohort,
      levels = c("ALSPAC", "UKB", "RS")
    )
  ) |>
  nest(Dat = -grp2) |>
  mutate(
    Plt = map(
      Dat,
      ~.x |>
        ggplot(aes(bmiq, estimate)) +
        geom_hline(yintercept = 0, lty = "11", linewidth = .2) +
        geom_linerange(
          aes(ymin = conf.low, ymax = conf.high, group = cohort),
          linewidth = 0.3, position = position_dodge(width = 0.25)
        ) +
        geom_point(
          aes(color = cohort, group = cohort),
          size = .5, position = position_dodge(width = 0.25)
        ) +
        facet_wrap(~sex) +
        ggh4x::facet_nested(
          grp1 + trait ~ sex, scales = "free_y",
          strip = ggh4x::strip_nested(size = "variable")
        ) +
        labs(
          x = "BMI percentiles",
          y = "Difference between\n75th - 25th slope percentiles",
          color = NULL
        ) +
        theme_bw()
    ) |>
      map2(
        grp2,
        ~ if (.y == 2) {
          .x + theme(axis.title.y = element_blank())
        } else {
          .x
        }
      )
  ) |>
  pull(Plt) |>
  patchwork::wrap_plots(nrow = 1) &
  theme(
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = -10, b = -10),
    legend.key.width = unit(1, "mm"),
    # Strips
    strip.text = element_text(size = 5, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.5, "mm"),
    # Axes
    axis.text.x = element_text(size = 4),
    axis.text.y = element_text(size = 2.5),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

outbmsdifs |>
  filter(trait %in% c("AndF_perc", "GynF_perc")) |>
  head()

options(repr.plot.width = 2, repr.plot.height = 3.5)
outbmsdifs |>
  filter(
    trait %in% c("GynF_perc", "AndF_perc"),
    cohort %in% c("RS", "UKB")
  ) |>
  pivot_wider(
    names_from = trait,
    values_from = c(estimate, conf.low, conf.high)
  ) |>
  mutate(grp1 = "Current BMI percentiles") |>
  ggplot(aes(estimate_AndF_perc, estimate_GynF_perc)) +
  geom_hline(yintercept = 0, lty = "11", linewidth = .2) +
  geom_vline(xintercept = 0, lty = "11", linewidth = .2) +
  # Fill area above the diagonal line
  geom_abline(slope = 1, intercept = 0, lty = "11", linewidth = .2) +
  geom_ribbon(
    stat = "function", fun = function(x) x,
    mapping = aes(ymin = after_stat(y), xmin = after_stat(x)),
    ymax = Inf,
    fill = "#41b9e1", alpha = 0.25
  ) +
  geom_ribbon(
    stat = "function", fun = function(x) x,
    mapping = aes(ymax = after_stat(y)),
    ymin = -Inf,
    fill = "#9a0740", alpha = 0.25
  ) +
  coord_equal() +
  geom_linerange(
    aes(xmin = conf.low_AndF_perc, xmax = conf.high_AndF_perc),
    linewidth = 0.3
  ) +
  geom_linerange(
    aes(ymin = conf.low_GynF_perc, ymax = conf.high_GynF_perc),
    linewidth = 0.3
  ) +
  geom_point(
    aes(fill = sex, shape = cohort)
  ) +
  scale_shape_manual(values = c(22, 23, 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  ggh4x::facet_nested(grp1 + factor(bmiq, levels = rev(unique(bmiq))) ~ .) +
  theme_bw() +
  labs(
    x = "∆AndF% due to higher BMI slope",
    y = "∆GynF% due to higher BMI slope",
    color = "Sex",
    shape = "Cohort"
  ) +
  theme(
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = -10, b = -10),
    legend.key.width = unit(1, "mm"),
    # Organize legends so they go vertically
    legend.box = "vertical",
    legend.spacing.y = unit(5, "mm"),
    # Strips
    strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
    panel.spacing = unit(0.25, "mm"),
    # Axes
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(linewidth = 0.2)
  )

addoutbms <- predtab |>

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
  filter(!trait %in% bmrks)
print(addoutbms)

options(repr.plot.width = 5, repr.plot.height = 3)
addoutbms |>
  filter(trait == "FI") |>
  ggplot(aes(bmi, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(
    aes(ymin = conf.low, ymax = conf.high)
  ) +
  geom_point(
    aes(color = cohort)
  ) +
  facet_wrap(~sex)

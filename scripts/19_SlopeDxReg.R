library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(survival)
library(mgcv)
library(patchwork)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

lmmtabrec <- read_rds(
  file.path(projfld, "data", "UKB", "processed_data", "lmmtabrec.rds")
)

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)

cvrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "cvrec.tsv"),
  show_col_types = FALSE
)

ancestrytab <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "ancestrytab.tsv"),
  col_types = "nf"
) |>
  mutate(ancestry = 1 * (ancestry == "EUR"))

outdat <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv"),
  show_col_types = FALSE
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
  # Preserve the existing survival adjustment, not a WHR outcome analysis.
  "WHR"
)

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

slopedxregdf <- bmrec |>
  # Joining biomarker and covariate data
  inner_join(cvrec, by = "eid") |>
  inner_join(ancestrytab, by = "eid") |>
  select(
    eid, sex, age, bmi,
    all_of(covars),
    all_of(bmrks)
  ) |>
  # Joining with survival data
  inner_join(outdat, by = "eid") |>
  # Nesting by sex and outcome
  nest(Dat = -c(sex, outcome)) |>
  # Joining with trajectory model data
  inner_join(lmmtabrec, by = "sex") |>
  # Joining with trajectory estimates
  transmute(
    sex, outcome,
    Dat = map2(indiv_est, Dat, inner_join, by = "eid")
  )

slopedxregdf <- slopedxregdf |>
  mutate(
    survcox_slopeonly = map(
      Dat,
      ~ coxph(
        Surv(time, event) ~ Slope,
        data = .x
      )
    )
  )

slopedxregdf <- slopedxregdf |>
  mutate(
    survcox_slopecat = map(
      Dat,
      ~ {
        xcat <- .x |>
          mutate(
            SlopeCat = ntile(Slope, 5) |>
              factor(labels = paste0("Q", 1:5)) |>
              relevel(ref = "Q3")
          )
        coxph(
          Surv(time, event) ~ SlopeCat,
          data = xcat
        )
      }
    )
  )

marginal_outcomes_plot <- slopedxregdf |>
  transmute(
    sex,
    outcome = factor(
      outcome,
      levels = c("CAD", "Stroke", "LiverD"),
      labels = c("CAD", "Stroke", "Liver disease")
    ),
    tidy_res = map(
      survcox_slopecat,
      ~ .x |>
        broom::tidy(conf.int = TRUE, exponentiate = TRUE) |>
        mutate(term = gsub("SlopeCat", "", term)) |>
        add_row(
          term = "Q3",
          estimate = 1, conf.low = 1, conf.high = 1
        )
    )
  ) |>
  unnest(tidy_res) |>
  ggplot(aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_point(
    aes(color = term == "Q3"),
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  facet_grid(sex ~ outcome) +
  theme_bw() +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "Slope quantile"
  )

saveRDS(
  marginal_outcomes_plot,
  file.path(projfld, "data", "UKB", "processed_data", "plot_marginal_outcomes.rds")
)

slopedxregdf <- slopedxregdf |>
  mutate(
    survcox_full = map(
      Dat,
      ~ coxph(
        Surv(time, event) ~ . - eid,
        data = .x
      )
    )
  )

slopedxregdf |>
  transmute(
    sex, outcome,
    tidy_res = map(
      survcox_full,
      ~ broom::tidy(.x, conf.int = TRUE, exponentiate = TRUE)
    )
  ) |>
  unnest(tidy_res) |>
  filter(term == "Slope")

smterms <- c(
  "s(age, bs = \"cr\")",   # Non-linear effect of age
  "s(bmi, bs = \"cr\")",   # Non-linear effect of BMI
  "ti(age, bmi)",        # Variation in the effect of BMI due to age
  "s(Slope, bs = \"cr\")", # Non-linear effect of the BMI slope
  "ti(bmi, Slope)",      # Variation in the effect of BMI due to slope
  "ti(Slope, age)"       # Variation in the effect of the slope due to age
)

slopedxregdf <- slopedxregdf |>
  mutate(
    survgam_fullsm = map(
      Dat,
      ~ gam(
        reformulate(
          response = "time",
          termlabels = c(
            smterms, covars, bmrks
          )
        ),
        data = .x,
        family = cox.ph(),
        weights = event
      )
    )
  )

slopedxregdf <- slopedxregdf |>
  mutate(
    sm_tb = map(
      survgam_fullsm,
      function(mod) {
        tibble(
          sm_pos = seq_along(smterms),
          sm_term = smterms
        ) |>
          mutate(
            sm_dat = map(sm_pos, \(x) mod$smooth[[x]]),
            sm_label = map_chr(sm_dat, \(x) x$label),
            sm_dim = map_int(sm_dat, \(x) x$dim),
            sm_labels = map2(
              sm_label, sm_dim,
              \(lab, dim) {
                if (dim > 1) {
                  paste0(lab, seq_len(dim))
                } else {
                  lab
                }
              }
            )
          ) |>
          select(sm_term, sm_labels) |>
          unnest(sm_labels) |>
          mutate(
            sm_pen = map_dbl(sm_labels, \(x) mod$sp[x])
          )
      }
    )
  )

slopedxregdf <- slopedxregdf |>
  mutate(
    survgam_noslope = map2(
      Dat, sm_tb,
      ~ {
        smsub <- smterms[1:3]
        smsubp <- filter(.y, sm_term %in% smsub)$sm_pen
        gam(
          reformulate(
            response = "time",
            termlabels = c(
              smsub,
              bmrks, covars
            )
          ),
          data = .x,
          family = cox.ph(),
          weights = event,
          sp = smsubp
        )
      }
    )
  )

outcome_slope_tests <- slopedxregdf |>
  transmute(
    sex, outcome,
    slp_test = map2_dbl(
      survgam_fullsm, survgam_noslope,
      ~ anova(.x, .y, test = "LRT")[2, 5]
    )
  ) |>
  mutate(padj = p.adjust(slp_test, method = "fdr"))
write_csv(
  outcome_slope_tests,
  file.path(projfld, "data", "UKB", "results", "outcome_slope_tests.csv")
)
saveRDS(
  slopedxregdf,
  file.path(projfld, "data", "UKB", "processed_data",
            "survival_models_current.rds")
)

slopedxregdf |>
  filter(sex == "Male", outcome == "Stroke") |>
  transmute(
    sex, outcome,
    slp_test = pmap(
      list(Dat, survgam_fullsm, survgam_noslope),
      function(dat, fullmod, noslopemod) {
        xt <- mutate(dat, time = 10)
        p1 <- -predict.gam(noslopemod, newdata = xt, type = "link")
        p2 <- -predict.gam(fullmod, newdata = xt, type = "link")
        concordance(
          Surv(time, event) ~ p1 + p2,
          data = dat
        )
      }
    )
  ) |>
  pluck("slp_test", 1) %>%
  {
    print(.)
    cntvec <- c(-1, 1)
    cdif <- cntvec %*% coef(.)
    cdifv <- sqrt(cntvec %*% vcov(.) %*% cntvec)
    ci <- qnorm(1 - 0.05 / 2) * cdifv
    lwr <- cdif - ci
    upr <- cdif + ci
    write_csv(
      tibble(
        sex = "Male", outcome = "Stroke",
        c_index_without_slope = coef(.)[1], c_index_with_slope = coef(.)[2],
        difference = as.numeric(cdif), conf_low = as.numeric(lwr),
        conf_high = as.numeric(upr)
      ),
      file.path(projfld, "data", "UKB", "results", "stroke_concordance.csv")
    )
    message(
      "\n",
      "C-indices", paste(round(coef(.), 3), collapse = ", "), "\n",
      "Concordance difference: ", round(cdif, 3), "\n",
      "CI: ", round(lwr, 3), " - ", round(upr, 3), "\n",
      "P-value:", round(pnorm(-abs(cdif / cdifv)) * 2, 3)
    )
  }

calibdat <- slopedxregdf |>
  filter(sex == "Male", outcome == "Stroke") |>
  mutate(
    across(
      c(survgam_fullsm, survgam_noslope),
      function(mod) {
        map2(
          mod, Dat,
          function(m, d) {
            xnew <- d
            xnew$time <- 10
            prob <- 1 - predict(m, newdata = xnew, type = "response")
            probcat <- ntile(prob, 10)
            srate <- map(
              1:10,
              ~{
                dcat <- d[which(probcat == .x), ]
                tryCatch(
                  {
                    sf <- survfit(Surv(time, event) ~ 1, data = dcat) |>
                      summary(times = 10)
                    tibble(
                      obs_rate = ifelse(length(sf$surv), 1 - sf$surv, 0),
                      lower = ifelse(length(sf$upper), 1 - sf$upper, 0),
                      upper = ifelse(length(sf$lower), 1 - sf$lower, 0)
                    )
                  },
                  error = function(e) {
                    tibble(obs_rate = 0, lower = 0, upper = 0)
                  }
                ) |>
                  mutate(
                    medprob = median(prob[which(probcat == .x)]),
                    botprob = quantile(prob[which(probcat == .x)], 0.025),
                    topprob = quantile(prob[which(probcat == .x)], 0.975)
                  )
              }
            ) |>
              bind_rows()
          }
        )
      }
    ),
    res = map2(
      survgam_fullsm, survgam_noslope,
      ~bind_rows(
        mutate(.x, mod = "full"),
        mutate(.y, mod = "noslope")
      )
    )
  ) |>
  select(res) |>
  unnest(res)

stroke_calibration_plot <- calibdat |>
  ggplot(aes(medprob, obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_line(aes(group = mod, color = mod)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = mod)) +
  geom_point(aes(color = mod)) +
  labs(
    x = "Predicted risk",
    y = "Observed event rate",
    color = "Model"
  ) +
  scale_color_discrete(
    labels = c("full" = "With slope", "noslope" = "Without slope")
  ) +
  theme_bw()

ggsave(
  file.path(projfld, "plots", "stroke_calibration.png"),
  plot = stroke_calibration_plot, width = 5, height = 4,
  units = "in", dpi = 300, bg = "white"
)
write_csv(
  calibdat,
  file.path(projfld, "data", "UKB", "results", "table_stroke_calibration.csv")
)

compquint_isol <- slopedxregdf |>
  transmute(
    sex, outcome,
    tidy_res = map2(
      Dat, survgam_fullsm,
      function(x, mod) {
        xcat <- x |>
          mutate(SlopeCat = ntile(Slope, 5)) |>
          summarise(meanSl = mean(Slope, na.rm = TRUE), .by = SlopeCat) |>
          mutate(SlRef = meanSl[SlopeCat == 3]) |>
          filter(SlopeCat != 3)
        meanx <- x |>
          summarise(
            across(
              -c(eid, Slope, time, event),
              ~mean(.x, na.rm = TRUE)
            )
          ) |>
          mutate(eid = 0, Slope = 0, time = 10)
        xcat |>
          mutate(
            compres = map2(
              meanSl, SlRef,
              function(slp, slpref) {
                marginaleffects::comparisons(
                  mod,
                  newdata = meanx,
                  variable = list(
                    Slope = data.frame(high = slp, low = slpref)
                  ),
                  type = "link"
                ) |>
                  mutate(
                    across(c(estimate, conf.low, conf.high), exp)
                  )
              }
            )
          ) |>
          unnest(compres) |>
          add_row(SlopeCat = 3, estimate = 1, conf.low = 1, conf.high = 1)
      }
    )
  ) |>
  unnest(tidy_res)

compare_subjects <- function(model, alt, ref) {
  # Get design matrix for both subjects
  x1 <- predict(model, newdata = ref, type = "lpmatrix")
  x2 <- predict(model, newdata = alt, type = "lpmatrix")
  # Contrast
  contr <- x2[1, ] - x1[1, ]
  # Covariance matrix
  vcmat <- insight::get_varcov(model)
  # Difference in predictions
  estimate <- contr %*% coef(model)
  # Calculate variance of difference
  var_diff <- contr %*% vcmat %*% contr
  se <- sqrt(var_diff[1, 1])
  df <- insight::get_df(model)
  t_crit <- qt(1 - .05 / 2, df = df)
  ci <- t_crit * se
  data.frame(
    estimate = as.numeric(estimate),
    conf.low = as.numeric(estimate) - ci,
    conf.high = as.numeric(estimate) + ci
  )
}

compquint_comb <- slopedxregdf |>
  transmute(
    sex, outcome,
    tidy_res = map2(
      Dat, survgam_fullsm,
      function(x, mod) {
        xcat <- x |>
          mutate(SlopeCat = ntile(Slope, 5)) |>
          summarise(
            across(-c(eid, time, event), ~mean(.x, na.rm = TRUE)),
            .by = SlopeCat
          ) |>
          mutate(eid = 0, time = 10)
        map(
          (1:5)[-3],
          function(slc) {
            compare_subjects(
              mod,
              alt = filter(xcat, SlopeCat == slc),
              ref = filter(xcat, SlopeCat == 3)
            ) |>
              mutate(SlopeCat = slc)
          }
        ) |>
          bind_rows() |>
          mutate(across(c(estimate, conf.low, conf.high), exp)) |>
          add_row(SlopeCat = 3, estimate = 1, conf.low = 1, conf.high = 1)
      }
    )
  ) |>
  unnest(tidy_res)

compquint_strat <- slopedxregdf |>
  transmute(
    sex, outcome,
    tidy_res = map2(
      Dat, survgam_fullsm,
      function(x, mod) {
        xmean <- x |>
          summarise(
            across(
              -c(eid, Slope, time, event, age, bmi),
              ~mean(.x, na.rm = TRUE)
            )
          ) |>
          mutate(eid = 0, Slope = 0, time = 10)
        xcat <- x |>
          mutate(SlopeCat = ntile(Slope, 5)) |>
          summarise(
            across(Slope, ~mean(.x, na.rm = TRUE)),
            .by = SlopeCat
          )
        xgrid <- expand_grid(
          age = c(45, 55, 65),
          bmi = c(22.5, 27.5, 32.5, 37.5),
          SlopeCat = 1:5
        )
        xgrid |>
          nest(
            predx = -SlopeCat
          ) |>
          mutate(
            predx = pmap(
              list(predx, SlopeCat),
              function(d, slc) {
                if (slc == 3) {
                  return(
                    mutate(
                      d,
                      estimate = 1, conf.low = 1, conf.high = 1
                    )
                  )
                }
                d <- cross_join(d, xmean)
                dcomp <- d |>
                  mutate(
                    high = filter(xcat, SlopeCat == slc)$Slope,
                    low = filter(xcat, SlopeCat == 3)$Slope
                  ) |>
                  select(high, low)
                marginaleffects::comparisons(
                  mod,
                  newdata = d,
                  variables = list(
                    Slope = dcomp
                  ),
                  type = "link"
                ) |>
                  transmute(
                    age, bmi,
                    contrast = paste(dcomp$high, " - ", dcomp$low),
                    across(c(estimate, conf.low, conf.high), exp)
                  )
              }
            )
          ) |>
          unnest(predx)
      }
    )
  ) |>
  unnest(tidy_res)

compslope <- slopedxregdf |>
  transmute(
    sex, outcome,
    tidy_res = map2(
      Dat, survgam_fullsm,
      function(x, mod) {
        xgrid <- x |>
          summarise(
            across(
              -c(eid, time, event, bmi, age),
              ~mean(.x, na.rm = TRUE)
            )
          ) |>
          mutate(eid = 0, time = 10) |>
          expand_grid(
            bmi = c(22.5, 27.5, 32.5, 37.5),
            age = c(45, 55, 65)
          )
        marginaleffects::comparisons(
          mod,
          newdata = xgrid,
          type = "link",
          variable = list(
            Slope = c(0.5, -0.5)
          )
        )
      }
    )
  ) |>
  unnest(tidy_res) |>
  mutate(across(c(estimate, conf.low, conf.high), exp))

adjusted_outcomes_plot <- compslope |>
  ggplot(aes(x = bmi, y = estimate)) +
  geom_hline(yintercept = 1, linetype = "11") +
  geom_linerange(
    aes(ymin = conf.low, ymax = pmin(conf.high, 6))
  ) +
  geom_point() +
  geom_point(
    aes(y = ifelse(conf.high > 6, 6, NA)), shape = 17
  ) +
  ggh4x::facet_nested(outcome ~ sex + age, scales = "free_y") +
  labs(
    x = "BMI",
    y = "Hazard Ratio (95% CI)"
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

saveRDS(
  adjusted_outcomes_plot,
  file.path(
    projfld, "data", "UKB", "processed_data",
    "plot_adjusted_outcomes.rds"
  )
)
write_csv(
  compquint_isol,
  file.path(
    projfld, "data", "UKB", "results",
    "table_incident_outcomes_marginal_quintiles.csv"
  )
)
write_csv(
  compslope,
  file.path(
    projfld, "data", "UKB", "results",
    "table_incident_outcomes_adjusted_contrasts.csv"
  )
)

dcadat <- slopedxregdf |>
  filter(sex == "Male", outcome == "Stroke") |>
  transmute(
    sex, outcome,
    dca_res = pmap(
      list(Dat, survgam_fullsm, survgam_noslope),
      function(x, modf, modn) {
        xnew <- x
        xnew$time <- 10
        x$p_modf <- 1 - predict(modf, newdata = xnew, type = "response")
        x$p_modn <- 1 - predict(modn, newdata = xnew, type = "response")
        x <- x[complete.cases(x$p_modf, x$p_modn), ]
        sv <- survfit(Surv(time, event) ~ 1, data = x)
        rt <- 1 - summary(sv, times = 10)$surv
        expand_grid(
          threshold = seq(0, 0.1, by = 0.001),
          label = c("p_modf", "p_modn", "none", "all")
        ) |>
          mutate(
            res = map2(
              threshold, label,
              function(thresh, lab) {
                if (lab == "none") {
                  return(tibble(rt = rt, ppos = 0, rt_in_pos = 0))
                }
                if (lab == "all") {
                  return(tibble(rt = rt, ppos = 1, rt_in_pos = rt))
                }
                pos <- which(x[[lab]] >= thresh)
                if (!length(pos)) {
                  return(tibble(rt = rt, ppos = 0, rt_in_pos = 0))
                }
                xpos <- x[pos, ]
                ppos <- nrow(xpos) / nrow(x)
                svpos <- tryCatch(
                  survfit(Surv(time, event) ~ 1, data = xpos),
                  error = function(e) NULL
                )
                rt_in_pos <- tryCatch(
                  1 - summary(svpos, times = 10)$surv,
                  error = function(e) 0
                )
                tibble(
                  rt = rt, ppos = ppos, rt_in_pos = rt_in_pos
                )
              }
            )
          ) |>
          unnest(res)
      }
    )
  )

dca_results <- dcadat |>
  select(sex, outcome, dca_res) |>
  unnest(dca_res) |>
  mutate(
    tpr = ppos * rt_in_pos,
    fpr = ppos * (1 - rt_in_pos),
    wt = threshold / (1 - threshold),
    nb_raw = tpr - (wt * fpr),
    nb = 100000 * nb_raw,
    nb_all = rt - (wt * (1 - rt)),
    nia = 100000 * ((nb_raw - nb_all) / wt)
  )

dca_nb_plot <- dca_results |>
  ggplot(aes(x = threshold, y = nb)) +
  coord_cartesian(ylim = c(-1e-3, NA)) +
  geom_smooth(
    aes(color = label),
    linewidth = .5,
    se = FALSE
  ) +
  scale_color_discrete(
    labels = c(
      "p_modf" = "With slope",
      "p_modn" = "Without slope",
      "all" = "Treat all",
      "none" = "Treat none"
    )
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Threshold probability",
    y = "Net benefit per 100k",
    color = "Model"
  ) +
  theme_bw()

dca_nia_plot <- dca_results |>
  ggplot(aes(x = threshold, y = nia, color = label)) +
  coord_cartesian(ylim = c(-1e-3, NA)) +
  geom_line() +
  scale_color_discrete(
    labels = c(
      "p_modf" = "With slope",
      "p_modn" = "Without slope",
      "all" = "Treat all",
      "none" = "Treat none"
    )
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Threshold probability",
    y = "Net interventions avoided per 100k",
    color = "Model"
  ) +
  theme_bw()

stroke_decision_plot <- patchwork::wrap_plots(
  dca_nb_plot, dca_nia_plot, nrow = 1
) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(projfld, "plots", "stroke_decision_curves.png"),
  plot = stroke_decision_plot, width = 10, height = 3,
  units = "in", dpi = 300, bg = "white"
)
write_csv(
  dca_results,
  file.path(
    projfld, "data", "UKB", "results",
    "table_stroke_decision_curves.csv"
  )
)

prophazfx <- function(model) {
  # Score residuals
  rs <- residuals(model, "score")
  # Term indices
  term_list <- attr(rs, "term")
  lapply(
    term_list,
    function(term) {
      cbind(
        time = model$y,
        rs[, term]
      ) |>
        data.frame() |>
        tidyr::pivot_longer(
          -time,
          names_to = "basis",
          values_to = "residual"
        )
    }
  ) |>
    dplyr::bind_rows(.id = "term")
}

phdat <- slopedxregdf |>
  filter(sex == "Male", outcome == "Stroke") |>
  pluck("survgam_fullsm", 1) |>
  prophazfx()

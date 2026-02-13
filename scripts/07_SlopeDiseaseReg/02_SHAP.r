library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(ggplot2)
library(survival)
library(randomForestSRC)
library(survex)

options(repr.plot.res = 300)

projfld <- "/ludc/Home/daniel_c/projects/DVA/LongitudinalBMI_SOPHIA"

bmrec <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "bmrec.tsv"),
  show_col_types = FALSE
)
head(bmrec)

outdat <- read_tsv(
  file.path(projfld, "data", "UKB", "processed_data", "outdat.tsv"),
  show_col_types = FALSE
)

survdat <- outdat |>
  filter(outcome == "CAD") |>
  select(-outcome) |>
  inner_join(bmrec, by = "eid") |>
  transmute(
    event, time,
    WHR, CRP, GLU, SBP,
    HDL,
    nonHDL = LDL + (TG / 2.2)
  ) |>
  drop_na() |>
  data.frame()
head(survdat)

coxmod <- coxph(
  Surv(time, event) ~ .,
  data = survdat,
  x = TRUE, model = TRUE
)

summary(coxmod)

rfmod <- rfsrc.fast(
  Surv(time, event) ~ .,
  data = survdat,
  ntree = 100,
  importance = "none",
  ntime = 20,
  perf.type = "none",
  forest = TRUE,
  save.memory = TRUE
)

preds <- randomForestSRC:::predict.rfsrc(
  rfmod, survdat,
  importance = "none",
  perf.type = "none",
)$survival[, 20]

concordance(Surv(survdat$time, survdat$event) ~ preds)

bg_ind <- sample(nrow(survdat), size = 1000, replace = FALSE)

expd <- survdat[bg_ind, ]

dx <- select(expd, -c(time, event))
ys <- with(expd, Surv(time, event))

explainer <- explain(
  model = rfmod,
  times = seq(0, 10, by = 2),
  data = dx,
  y = ys
)

x_ev <- survdat |>
  filter(event == 1, time > 9)

x_ev_predx <- mutate(x_ev, time = 10, event = 0)
x_ev$pred <- 1 - predict(coxmod, newdata = x_ev_predx, type = "survival")
x_ev <- x_ev |>
  filter(pred >= 0.10 & pred <= 0.15)

x_evd <- select(x_ev, -c(event, time, pred))

nrow(x_evd)

shapvals <- kernelshap::kernelshap(
  coxmod,
  X = x_evd,
  bg_X = dx,
  exact = TRUE,
  pred_fun = survival:::predict.coxph
)

shapv <- shapvals$S |>
  data.frame() |>
  mutate(id = row_number()) |>
  pivot_longer(-id, names_to = "variable", values_to = "shap_value")

options(repr.plot.width = 4, repr.plot.height = 4)
shapv |>
  ggplot(aes(shap_value, variable)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Change in logHR", y = NULL)

exemplars <- shapv |>
  group_by(variable) |>
  summarise(
    max_shap = max(shap_value),
    max_shapidx = first(id[shap_value == max_shap]),
  ) |>
  pivot_longer(
    -variable,
    names_to = c("label", ".value"),
    names_sep = "_"
  ) |>
  mutate(
    dat = map(
      shapidx,
      ~x_ev[.x, ]
    )
  ) |>
  unnest(dat)
exemplars

shapdat <- predict_parts(
  explainer,
  new_observation = exemplars[, -c(1:6)],
  type = "survshap",
  output_type = "survival",
  aggregation_method = NULL,
)

extds <- exemplars |>
  mutate(
    id = row_number(),
    tdshap = map(
      id,
      ~shapdat |>
        pluck("result", .x) |>
        tibble::rownames_to_column(var = "time") |>
        mutate(time = as.numeric(gsub("t=", "", time)))
    )
  )
print(extds)

toplot <- extds |>
  rename(time_ev = time) |>
  mutate(
    tdshap = map(tdshap, ~rename_with(.x, ~paste0(.x, "_shap"), -time))
  ) |>
  unnest(tdshap) |>
  pivot_longer(
    ends_with("_shap"),
    names_to = "feature", values_to = "shapval"
  ) |>
  mutate(feature = gsub("_shap$", "", feature))
head(toplot)

options(repr.plot.width = 7, repr.plot.height = 5)
toplot |>
  ggplot(aes(x = time, y = shapval * 100)) +
  geom_line(aes(group = feature, color = feature)) +
  guides(color = guide_legend(nrow = 1)) +
  geom_text(
    data = . %>%
      group_by(shapidx) %>%
      summarise(
        label = paste0(
          "- WHR = ", unique(round(WHR, 1)),
          "\n- CRP = ", unique(round(CRP, 1)), " mg/L",
          "\n- GLU = ", unique(round(GLU, 1)), " mmol/L",
          "\n- SBP = ", unique(round(SBP, 1)), " mmHg",
          "\n- HDL = ", unique(round(HDL, 1)), " mmol/L",
          "\n- nonHDL = ", unique(round(nonHDL, 1)), " mmol/L"
        )
      ),
    aes(x = 0.5, y = -5, label = label),
    hjust = 0, size = 1.75
  ) +
  facet_wrap(~shapidx) +
  labs(
    x = "Time (years)",
    y = "% change in predicted\nsurvival probability",
    color = "Feature"
  ) +
  theme_bw() +
  theme(legend.position = "top")

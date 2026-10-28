# Recent BMI trajectories and metabolic profiles at the same BMI

Analysis code for the paper on recent BMI slope and metabolic profiles at
the same BMI across three life stages: ALSPAC (young adults), UK Biobank
(midlife) and the Rotterdam Study (older adults).

Individual-level data are not shared here. UK Biobank data are available to
approved researchers (this work used application 57232); ALSPAC and Rotterdam
analyses were run by each cohort's analysts following
`scripts/21_Guideline.ipynb`, and only model results were shared.

## Contents

- `scripts/`: analysis scripts, numbered by stage. `src/run_pipeline.sh`
  lists the order in which they are run. Each script has a header describing
  its purpose, inputs and outputs.
- `src/`: helper functions and the pipeline runner.

## Analysis outline

1. Data preparation and 5-year BMI slopes from primary-care weight records
   (scripts 01–04, 11).
2. Generalized additive models of biomarkers, DXA fat and liver fat on age,
   BMI and slope (05, 12, 13, 22); subsampling and sensitivity analyses
   (06–10, 14–17, 32).
3. Within-cohort contrasts between the lower and upper quartile of slope at
   the same BMI, evaluated from the fitted models (36–38); waist-to-hip ratio
   follow-up for liver fat (35, 46–48).
4. Gain, stable and loss contrasts in UK Biobank and Rotterdam (39, 40).
5. Lagged biomarkers, incident outcomes and the Diabetes Prevention Program
   comparison (18–20).
6. Cross-cohort summaries, figures, supplementary tables and documents
   (23–31, 41–45); table and document assembly is in `paper/`.

## Software

R 4.5.2 with mgcv 1.9-4, dplyr, tidyr, purrr, readr, lme4, survival,
ggplot2 and patchwork. Python 3 (matplotlib, pandas) for diagrams.

## License

GPL-3.0; see `LICENSE`.

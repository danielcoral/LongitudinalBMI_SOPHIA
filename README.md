# Recent BMI trajectories and metabolic profiles at the same BMI

Analysis code for the paper on recent BMI slope and metabolic profiles at
the same BMI across three life stages: ALSPAC (young adults), UK Biobank
(midlife) and the Rotterdam Study (older adults).

Individual-level data are not shared here. UK Biobank data are available to
approved researchers (this work used application 57232); ALSPAC and Rotterdam
analyses were run by each cohort's analysts following
`scripts/21_Guideline.ipynb`, and only model results were shared.

## Contents

- `scripts/`: numbered analysis stages; see `scripts/README.md` for setup,
  stage groups, workflow commands and HTCondor execution.
- `src/`: reusable analysis and figure helpers.
- `config/config.yaml`: shared paths, cohort result locations and each
  stage's inputs, outputs and resource settings.
- `Snakefile`: dependency-aware analysis entry point.
- `data/local/{raw,processed}/`: restricted UK Biobank inputs and derived
  participant-level datasets, saved models and plotting objects.
- `data/external/processed/`: returned ALSPAC and Rotterdam model results.
- `results/{ukb,cross_cohort}/` and `plots/`: generated analysis exports.
- `paper/`: manuscript, bibliography, asset maps and document tools; see
  `paper/README.md` before building publication assets.

Data, results and publication assets are private and excluded from Git.
Keep existing received data unchanged; model outputs can contain sensitive
information even when they are not rectangular participant tables.

## Analysis outline

1. Prepare data and 5-year BMI slopes from primary-care weight records
   (01-04, 11).
2. Fit generalized additive models of biomarkers, DXA fat and liver fat on
   age, BMI and slope (05, 12, 13, 22); subsampling and sensitivity analyses
   (06-10, 14-17, 32).
3. Estimate within-cohort lower-versus-upper slope-quartile contrasts at the
   same BMI (36-38); evaluate liver-fat/waist-to-hip-ratio follow-ups
   (35, 46-48).
4. Compare gain, stable and loss contrasts in UK Biobank and Rotterdam
   (33, 34, 39, 40).
5. Evaluate lagged biomarkers, incident outcomes and the Diabetes Prevention
   Program comparison (18-20).
6. Assemble cross-cohort summaries, figures and manuscript numbers
   (23-29, 31, 41-45); build publication exports separately in `paper/`.

Script numbers identify stages; the Snakefile resolves execution order
from their saved-file dependencies.

## Software and execution

`environment.yml` records the Conda environment and direct requirements,
including the additional CRAN installation command. The existing analysis
uses R 4.5.2 and mgcv 1.9-4. R figures use ggplot2; two existing Python
diagrams use matplotlib.

With the project environment active, inspect the workflow from this root:

```sh
snakemake --dry-run --cores 1
```

Full analyses run through HTCondor; setup, submission and verification
commands are in `scripts/README.md`. Publication builds are explicit
opt-in targets and are excluded from the default analysis run.

## License

GPL-3.0; see `LICENSE`.

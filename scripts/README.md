# Analysis workflow

Run from the project root. `config/config.yaml` owns shared paths and each
stage's inputs, outputs and resource requests; `src/project_setup.R` loads
those paths for R scripts. `Snakefile` resolves the execution order.

## Environment

The existing RLang environment is retained. `environment.yml` is its single
dependency manifest; no environment was changed during template mapping.
Activate the environment using your own Miniforge installation:

```sh
conda activate RLang
```

On another machine, create it with `conda env create --file environment.yml`,
then run the additional CRAN installation command recorded in that file.
Record package changes in the manifest before installing them. Use
`conda env update --file environment.yml` only for intentional changes;
the existing environment may also support other projects.

Select this environment with VS Code's **Python: Select Interpreter**.
All commands below assume it is active.

## Stages and inputs

- 01: optional GP weight extraction from `data/local/raw/gp_clinical.csv`.
  The existing `OnlyWeights_gp.tsv` is an input to the default workflow;
  extraction is not automatically repeated.
- 02-05, 11-13: data preparation, trajectories and primary/imaging/liver-fat
  models. Participant-level data and fitted objects stay in
  `data/local/processed/`.
- 06-10, 14-17, 32: subsampling and sensitivity analyses; recruitment and
  imaging bootstrap files have explicit dependencies.
- 18-20, 24: lagged biomarkers, incident outcomes, DPP and outcome figure.
- 21: standalone collaborator guideline notebook, retained unchanged;
  `handoffs/Guideline.html` is the existing rendered guideline. Neither is
  executed by Snakemake.
- 22-23, 27-29: UKB cohort export, cross-cohort summaries, matched contrasts
  and manuscript-number exports. External ALSPAC/Rotterdam model returns
  are inputs configured under `cohort_results`.
- 30: optional alternative reporting from saved fits; excluded from the DAG
  because it overlaps the outputs of 05/12/13.
- 31: final composites from saved panels, after the producing model and
  subsampling stages.
- 33-41: stable-reference/asymmetry, within-cohort/support and liver-fat
  contrasts, WHR interaction and results digest.
- 42-45: main figures from saved contrasts and summaries.
- 46-48: separate WHR models, their summaries and liver-fat WHR adjustment.
- 25-26: optional Python diagrams with fixed study wording/counts. Review
  those annotations against current results before using regenerated images.

Numbering is retained for continuity. In particular, 36 reads outputs from
27 and 33; 37 reads summaries from 29; 31 uses saved output from 17.
Stage 17 combines raw-liver subsamples from 14 with log-liver subsamples
from 16.
The workflow declares those dependencies rather than sorting filenames.
Some original model scripts also write preliminary figures later replaced
by 31; the final composite belongs to 31 in the DAG.

## Inspect and run

Inspection does not run analyses:

```sh
snakemake --dry-run --cores 1
snakemake --list-rules
```

For an authorized full run, the established execution mode is HTCondor:

```sh
condor_submit subfiles/pipeline.sub
condor_q JOB_ID
condor_history JOB_ID -limit 1
```

The job runs Snakemake on the shared filesystem with four CPUs and 64 GB
memory. Per-stage requests are conservative starting settings; assess real
usage before increasing concurrency. Snakemake skips current outputs.
It tracks source/configuration changes, so the first run after migration
can require refits even though saved outputs were preserved.

Scheduler logs use `logfiles/pipeline.<cluster>.*`; stage logs use
`logfiles/stage_<number>.log`. Do not launch overlapping full/individual
jobs against the same outputs. Check exit status, stage logs and expected
artifacts before marking a run complete. Jobs were not submitted during
mapping.

`bash src/run_pipeline.sh --dry-run` retains the old entry point for
inspection. For an individual long stage, prepare a copy of
`subfiles/SUBTEMPLATE.sub` with absolute project/Conda paths, the environment
name from `environment.yml`, resources, logs and `initialdir` set to this
root. The three existing bootstrap submissions remain available.

## Optional exports and checks

These are explicit Snakemake targets, outside the default analysis:
`diagrams`, `main_figures`, `supplementary_tables`, `manuscript`,
`supplementary_notes` and `supplementary_figures`. Inspect them with
`--dry-run`. Actual document
builds overwrite their outputs, so read `paper/README.md` first.
No missing legacy `paper/assemble_supplementary_tables.R` or
`paper/rebuild_manuscript_documents.sh` is called by the new workflow.

Focused checks use synthetic inputs and do not read participant data:

```sh
Rscript --vanilla tests/test_project_mapping.R
Rscript --vanilla paper/scripts/test_builders.R
python paper/scripts/test_suppfigs.py
python paper/scripts/test_extract_review.py
```

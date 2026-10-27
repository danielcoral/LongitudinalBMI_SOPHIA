#!/usr/bin/env bash
set -euo pipefail

project="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
conda="/ludc/Home/daniel_c/miniforge3/bin/conda"
run_log="$project/logfiles/pipeline_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "$project/logfiles"
cd "$project"

scripts=(
  02_DataPrep.R
  03_Descriptive.R
  04_Traj_rec.R
  # Shared biomarkers and DXA outcomes across cohorts.
  05_GAMRec.R
  11_Traj_Img2i.R
  12_GAMImg2i.R
  22_FollowGuideUKB.R
  23_ReplicationRes.R
  37_WithinCohortContrasts.R
  38_SlopeSupportSummary.R
  # Liver fat, then imaging-visit WHR as an external-measure follow-up.
  13_GAMlogLF.R
  39_LiverFatQuartileContrasts.R
  36_LiverFatWhrInteraction.R
  # Stable-reference analyses use shared traits and liver fat, not WHR.
  34_GainVsStableContrasts.R
  35_SummarizeSlopeAsymmetry.R
  40_StableContrastsCommonGrid.R
  41_StableConsistencyUkbRs.R
  42_ResultsDigest.R
  # Sensitivity/subsampling: main biomarkers, DXA and liver fat only.
  06_GAMRec_Boot.R
  07_GAMRec_BootRes.R
  14_GAMImg2i_Boot.R
  15_GAMImg2i_BootRes.R
  16_GAMlogLF_Boot.R
  17_GAMlogLF_BootRes.R
  08_NoDxRec.R
  09_StrictFollowupRec.R
  10_SecondVisitRec.R
  27_MatchedCohortContrasts.R
  28_PlotMatchedCohorts.R
  33_SummarizeSensitivityAgreement.R
  18_SlopeBMRegDelay.R
  19_SlopeDxReg.R
  20_RCTComp.R
  24_AssembleOutcomesSummary.R
  # Reporting and document assembly.
  29_ExportManuscriptNumbers.R
  31_AssembleSavedFigures.R
  43_PlotDesignSlopes.R
  44_PlotWithinCohortContrasts.R
  45_PlotLiverFat.R
  46_PlotLiverFatWhrProfiles.R
  47_PlotStableReference.R
  32_AssembleSupplementaryTables.R
)

printf 'Pipeline started: %s\n' "$(date --iso-8601=seconds)" | tee -a "$run_log"
printf 'Excluded by request: 01_ExtractGPWeight.sh, 21_Guideline.ipynb\n' | tee -a "$run_log"
for script in "${scripts[@]}"; do
  printf '\n[%s] START scripts/%s\n' "$(date --iso-8601=seconds)" "$script" | tee -a "$run_log"
  "$conda" run --name RLang --no-capture-output Rscript --vanilla \
    "scripts/$script" >>"$run_log" 2>&1
  printf '[%s] DONE  scripts/%s\n' "$(date --iso-8601=seconds)" "$script" | tee -a "$run_log"
done
"$conda" run --name RLang --no-capture-output bash \
  paper/rebuild_manuscript_documents.sh >>"$run_log" 2>&1
printf '\nPipeline completed: %s\n' "$(date --iso-8601=seconds)" | tee -a "$run_log"

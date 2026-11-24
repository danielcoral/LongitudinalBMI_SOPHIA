#!/usr/bin/env bash
# Extract GP weight fields. Input: gp_clinical.csv. Output: OnlyWeights_gp.tsv.
set -euo pipefail
cd "$(dirname "$0")/.."
raw_dir=$(python - <<'PY'
import yaml
with open("config/config.yaml") as stream:
    print(yaml.safe_load(stream)["paths"]["ukb_raw"])
PY
)
temporary=$(mktemp "$raw_dir/OnlyWeights_gp.XXXXXX")
trap 'rm -f -- "$temporary"' EXIT
awk -F, 'BEGIN {OFS="\t"}
  NR == 1 || $5 == "22A.." || $6 == "22A.." {
    print $1, $3, $4, $7, $8, $9
  }' "$raw_dir/gp_clinical.csv" > "$temporary"
mv -- "$temporary" "$raw_dir/OnlyWeights_gp.tsv"

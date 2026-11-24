#!/usr/bin/env bash
# Compatibility entry point; the root Snakefile owns workflow dependencies.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
conda="${CONDA_EXE:-$HOME/miniforge3/bin/conda}"
environment=$(sed -n 's/^name: *//p' environment.yml)
exec "$conda" run --name "$environment" --no-capture-output \
  snakemake --cores "${SNAKEMAKE_CORES:-1}" "$@" \
  --resources mem_mb=64000

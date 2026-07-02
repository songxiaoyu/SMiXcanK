#!/usr/bin/env bash

# Run the HEARTFAIL PWAS workflow from liftover through plotting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"

export HEARTFAIL_N_CASES="${HEARTFAIL_N_CASES:-1405}"
export HEARTFAIL_N_CONTROLS="${HEARTFAIL_N_CONTROLS:-359789}"
export HEARTFAIL_ASSOC_REG_MODE="${HEARTFAIL_ASSOC_REG_MODE:-fixed}"
export HEARTFAIL_ASSOC_REG_SCALE="${HEARTFAIL_ASSOC_REG_SCALE:-0.1}"
export HEARTFAIL_SKIP_EXISTING="${HEARTFAIL_SKIP_EXISTING:-true}"
export HEARTFAIL_PARALLEL_JOBS="${HEARTFAIL_PARALLEL_JOBS:-4}"
export PLINK2_BIN="${PLINK2_BIN:-/opt/anaconda3/envs/plink2/bin/plink2}"

cd "${REPO_DIR}"

echo "Step 1: liftover HEARTFAIL GWAS"
python3 Analysis_code/HEARTFAIL/1_liftover_heartfail_gwas.py

echo "Step 2: prepare HEARTFAIL PWAS input"
"${RSCRIPT_BIN}" Analysis_code/HEARTFAIL/2_HEARTFAIL_prepare_data_pwas.R

echo "Step 3: build EUR LD reference"
bash Analysis_code/HEARTFAIL/3_1000Genome_keep_eur_plink_heartfail_pwas.sh

echo "Step 4: association"
bash Analysis_code/HEARTFAIL/run_HEARTFAIL_step4_parallel.sh

echo "Step 5: Primo"
"${RSCRIPT_BIN}" Analysis_code/HEARTFAIL/5_run_Primo_heartfail_pwas.R

echo "Step 6: plot"
"${RSCRIPT_BIN}" Analysis_code/HEARTFAIL/6_plot_heartfail_pwas.R

echo "Done."
echo "Result:"
echo "  /Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv"
echo "Figure dir:"
echo "  /Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Figure/heartfail_pwas"

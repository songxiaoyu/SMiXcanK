#!/usr/bin/env bash

# Run the HERMES HFpEF PWAS workflow from liftover through plotting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_hfpef_pwas/hermes_hfpef_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
FIGURE_DIR="${PAPER_DIR}/Figure/hermes_hfpef_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
RSCRIPT_BIN="/Library/Frameworks/R.framework/Resources/bin/Rscript"
PYTHON_BIN="/opt/anaconda3/bin/python"

cd "${REPO_DIR}"

echo "Step 1: liftover HERMES HFpEF GWAS"
"${PYTHON_BIN}" Analysis_code/HERMES_HFpEF/1_liftover_hermes_hfpef_gwas.py

echo "Step 2: prepare HERMES HFpEF PWAS input"
"${RSCRIPT_BIN}" Analysis_code/HERMES_HFpEF/2_HERMES_HFpEF_prepare_data_pwas.R

echo "Step 3: build EUR LD reference"
bash Analysis_code/HERMES_HFpEF/3_1000Genome_keep_eur_plink_hermes_hfpef_pwas.sh

echo "Step 4: association"
bash Analysis_code/HERMES_HFpEF/run_HERMES_HFpEF_step4_parallel.sh

echo "Step 5: Primo"
"${RSCRIPT_BIN}" Analysis_code/HERMES_HFpEF/5_run_Primo_hermes_hfpef_pwas.R

echo "Step 6: plot"
"${RSCRIPT_BIN}" Analysis_code/HERMES_HFpEF/6_plot_hermes_hfpef_pwas.R

echo "Done."
echo "Result: ${WORKSPACE_DIR}/hermes_hfpef_result/hermes_hfpef_result_pwas.csv"
echo "Annotated: ${WORKSPACE_DIR}/hermes_hfpef_result/hermes_hfpef_result_pwas_fixed_0p005_annotated.csv"
echo "Figure dir: ${FIGURE_DIR}"

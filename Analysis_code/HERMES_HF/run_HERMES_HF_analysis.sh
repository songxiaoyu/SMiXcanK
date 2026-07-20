#!/usr/bin/env bash

# Run the HERMES HF PWAS workflow from liftover through plotting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_hf_pwas/hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
FIGURE_DIR="${PAPER_DIR}/Figure/hermes_hf_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
RSCRIPT_BIN="/Library/Frameworks/R.framework/Resources/bin/Rscript"
PYTHON_BIN="/opt/anaconda3/bin/python"

cd "${REPO_DIR}"

echo "Step 1: liftover HERMES HF GWAS"
"${PYTHON_BIN}" Analysis_code/HERMES_HF/1_liftover_hermes_hf_gwas.py

echo "Step 2: prepare HERMES HF PWAS input"
"${RSCRIPT_BIN}" Analysis_code/HERMES_HF/2_HERMES_HF_prepare_data_pwas.R

echo "Step 3: build EUR LD reference"
bash Analysis_code/HERMES_HF/3_1000Genome_keep_eur_plink_hermes_hf_pwas.sh

echo "Step 4: association"
bash Analysis_code/HERMES_HF/run_HERMES_HF_step4_parallel.sh

echo "Step 5: Primo"
HERMES_HF_RESULT_TAG="fixed_0p200" "${RSCRIPT_BIN}" Analysis_code/HERMES_HF/5_run_Primo_hermes_hf_pwas.R

echo "Step 6: plot"
"${RSCRIPT_BIN}" Analysis_code/HERMES_HF/6_plot_hermes_hf_pwas.R

echo "Done."
echo "Result: ${WORKSPACE_DIR}/hermes_hf_result/hermes_hf_result_pwas_fixed_0p200.csv"
echo "Annotated: ${WORKSPACE_DIR}/hermes_hf_result/hermes_hf_result_pwas_fixed_0p200_annotated.csv"
echo "Table: ${WORKSPACE_DIR}/hermes_hf_result/hermes_hf_table_pwas_fixed_0p200.csv"
echo "Figure dir: ${FIGURE_DIR}"

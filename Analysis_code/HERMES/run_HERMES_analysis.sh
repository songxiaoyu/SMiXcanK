#!/usr/bin/env bash

# Run the current HERMES PWAS workflow with the restored moderate workspace.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
FIGURE_DIR="${PAPER_DIR}/Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
RSCRIPT_BIN="/Library/Frameworks/R.framework/Resources/bin/Rscript"

cd "${REPO_DIR}"

echo "Step 4: association"
bash Analysis_code/HERMES/run_HERMES_step4_parallel.sh

echo "Step 5: Primo"
"${RSCRIPT_BIN}" Analysis_code/HERMES/5_run_Primo_hermes_pwas.R

echo "Step 6: plot"
"${RSCRIPT_BIN}" Analysis_code/HERMES/6_plot_hermes_pwas.R

echo "Done."
echo "Result: ${WORKSPACE_DIR}/hermes_result/hermes_result_pwas.csv"
echo "Annotated: ${WORKSPACE_DIR}/hermes_result/hermes_result_pwas_fixed_0p005_annotated.csv"
echo "Figure dir: ${FIGURE_DIR}"

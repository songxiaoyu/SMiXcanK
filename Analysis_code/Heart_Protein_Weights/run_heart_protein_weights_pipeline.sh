#!/usr/bin/env bash

# Reproduce the current heart protein weight table from pruning through combine.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"

cd "${REPO_DIR}"

echo "Step 1: moderate genotype pruning"
PWAS_PRUNE_KB=100 \
PWAS_PRUNE_STEP=1 \
PWAS_PRUNE_R2=0.99 \
bash Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh

echo "Step 2: train chromosome-level weights"
bash Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh

echo "Step 3: combine chromosome-level weights"
"${RSCRIPT_BIN}" Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R

echo "Done."
echo "Output:"
echo "  /Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv"

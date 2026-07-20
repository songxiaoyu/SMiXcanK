#!/bin/bash

# Train heart protein weights chromosome by chromosome with the selected less-sparse setup:
#   genotype input: moderate pruning, 100kb window, r2 = 0.99
#   MiXcan alpha: 0.5
#   glmnet lambda: lambda.min
#
# This with-intercepts workflow writes new non-overwriting files, including:
#   weights_heart_protein_cardiomyocytes_other*_with_intercepts*.csv
#   intercepts_heart_protein_cardiomyocytes_other*_with_intercepts*.csv
#   model_terms_heart_protein_cardiomyocytes_other*_with_intercepts*.csv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"
PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
GENO_DIR="${PWAS_GENO_RAW_DIR:-${PAPER_DIR}/New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99}"
LOG_DIR="${PWAS_MODERATE_TRAIN_LOG_DIR:-/private/tmp/pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts_logs}"
CHR_LIST="${PWAS_CHR_LIST:-$(printf "%s," {1..22})}"
CHR_LIST="${CHR_LIST%,}"
OUTPUT_PREFIX="${PWAS_OUTPUT_PREFIX:-_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts}"
MAX_JOBS="${PWAS_MAX_JOBS:-4}"
SKIP_EXISTING="${PWAS_SKIP_EXISTING:-1}"

mkdir -p "${LOG_DIR}"
IFS=',' read -r -a CHRS <<< "${CHR_LIST}"

cd "${REPO_DIR}"

run_one_chr() {
  local chr="$1"
  local suffix="${OUTPUT_PREFIX}_chr${chr}"
  local results_dir="${PAPER_DIR}/Results/heart_protein_weights/training_model_weights"
  local weights_file="${results_dir}/weights_heart_protein_cardiomyocytes_other${suffix}.csv"
  local intercepts_file="${results_dir}/intercepts_heart_protein_cardiomyocytes_other${suffix}.csv"
  local model_terms_file="${results_dir}/model_terms_heart_protein_cardiomyocytes_other${suffix}.csv"

  if [ "${SKIP_EXISTING}" = "1" ] && [ -s "${weights_file}" ] && [ -s "${intercepts_file}" ] && [ -s "${model_terms_file}" ]; then
    echo "===== Skipping chr${chr}; complete outputs already exist ====="
    return 0
  fi

  echo "===== Training chr${chr} ====="
  PWAS_GENO_RAW_DIR="${GENO_DIR}" \
  PWAS_CHR_FILTER="${chr}" \
  PWAS_OUTPUT_SUFFIX="${suffix}" \
  PWAS_MIXCAN_ALPHA=0.5 \
  PWAS_LAMBDA_CHOICE=min \
  PWAS_WEIGHT_EPS=1e-8 \
  "${RSCRIPT_BIN}" Analysis_code/Heart_Protein_Weights/2b_train_heart_protein_weights_with_intercepts.R \
    > "${LOG_DIR}/train_chr${chr}.log" 2>&1
  echo "===== Finished chr${chr} ====="
}

echo "Training heart protein weights with moderate pruning, lambda.min, and intercept outputs"
echo "Genotype dir: ${GENO_DIR}"
echo "Chromosomes: ${CHR_LIST}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Logs: ${LOG_DIR}"
echo "Parallel jobs: ${MAX_JOBS}"
echo "Skip existing complete chromosome outputs: ${SKIP_EXISTING}"

pids=()
for chr in "${CHRS[@]}"; do
  while [ "$(jobs -rp | wc -l | tr -d '[:space:]')" -ge "${MAX_JOBS}" ]; do
    sleep 10
  done
  run_one_chr "${chr}" &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    status=1
  fi
done

if [ "${status}" -ne 0 ]; then
  echo "At least one chromosome failed. Check logs in ${LOG_DIR}" >&2
  exit "${status}"
fi

echo "All requested chromosomes finished."
echo "Combining outputs..."
PWAS_OUTPUT_PREFIX="${OUTPUT_PREFIX}" \
PWAS_CHR_LIST="${CHR_LIST}" \
"${RSCRIPT_BIN}" Analysis_code/Heart_Protein_Weights/3b_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights_with_intercepts.R

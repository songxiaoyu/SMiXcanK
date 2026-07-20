#!/usr/bin/env bash

# Run HERMES_HF HF PWAS step 4 for chr1-22, four chromosomes at a time.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_hf_pwas/hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
LOG_DIR="${WORKSPACE_DIR}/logs/step4_parallel"
RSCRIPT_BIN="/Library/Frameworks/R.framework/Resources/bin/Rscript"
JOBS=4
ASSOC_REG_SCALE="${HERMES_HF_ASSOC_REG_SCALE:-0.2}"
SCALE_TAG="$(printf "%.3f" "${ASSOC_REG_SCALE}" | tr "." "p")"
RESULT_TAG="fixed_${SCALE_TAG}"

mkdir -p "${WORKSPACE_DIR}/hermes_hf_result" "${LOG_DIR}"

# Each worker calls step 4 with one chromosome. Step 4 skips that chromosome
# when its result CSV already exists, so interrupted runs can be restarted.
run_chr() {
  local chr="$1"
  local log_file="${LOG_DIR}/chr${chr}.log"
  echo "Starting chr${chr}; log: ${log_file}"
  (
    cd "${REPO_DIR}"
    HERMES_HF_ASSOC_REG_SCALE="${ASSOC_REG_SCALE}" \
      "${RSCRIPT_BIN}" Analysis_code/HERMES_HF/4_HERMES_HF_run_analysis_pwas.R "${chr}"
  ) >"${log_file}" 2>&1
}

for chr in {1..22}; do
  while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "${JOBS}" ]; do
    sleep 5
  done
  run_chr "${chr}" &
done

status=0
for pid in $(jobs -pr); do
  if ! wait "${pid}"; then
    status=1
  fi
done

if [ "${status}" -ne 0 ]; then
  echo "At least one chromosome failed. Check logs in ${LOG_DIR}." >&2
  exit "${status}"
fi

echo "Combining chromosome results..."
(
  cd "${REPO_DIR}"
  HERMES_HF_RESULT_TAG="${RESULT_TAG}" "${RSCRIPT_BIN}" -e '
    result_tag <- Sys.getenv("HERMES_HF_RESULT_TAG", unset = "fixed_0p200")
    result_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/hermes_hf_pwas/hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_hf_result"
    files <- file.path(result_dir, sprintf("hermes_hf_chr%d_result_pwas_%s.csv", 1:22, result_tag))
    missing <- files[!file.exists(files)]
    if (length(missing)) {
      stop("Missing chromosome result files:\n", paste(missing, collapse = "\n"), call. = FALSE)
    }
    combined <- do.call(rbind, lapply(files, read.csv))
    out <- file.path(result_dir, sprintf("hermes_hf_result_pwas_%s.csv", result_tag))
    write.csv(combined, out, row.names = FALSE)
    cat("Wrote:", out, "\n")
    cat("Rows:", nrow(combined), "\n")
  '
)

echo "Done."

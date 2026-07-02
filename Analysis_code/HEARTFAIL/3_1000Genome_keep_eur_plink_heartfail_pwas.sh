#!/bin/bash

# Build 1000Genome EUR LD reference files for the HEARTFAIL PWAS SNP lists.
# Run after 2_HEARTFAIL_prepare_data_pwas.R.

set -euo pipefail

CONDA_SH="${CONDA_SH:-/opt/anaconda3/etc/profile.d/conda.sh}"
if [ -f "${CONDA_SH}" ]; then
  source "${CONDA_SH}"
fi

PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
DATA_DIR="${PAPER_DIR}/Data/plink_snplist_by_gene"
WORKSPACE_DIR="${HEARTFAIL_WORKSPACE_DIR:-${PAPER_DIR}/Results/heartfail_pwas/heartfail_workspace}"
GWAS_ID_DIR="${WORKSPACE_DIR}/heartfail_filtered_id"
EUR_SAMPLES="${PAPER_DIR}/Data/1000Genome/eur_ids.txt"

CHR_LIST="${HEARTFAIL_CHR_LIST:-$(echo {1..22})}"
PLINK2_ENV="${PWAS_PLINK2_ENV:-}"
PLINK2_BIN="${PLINK2_BIN:-plink2}"

for chr in ${CHR_LIST}; do
  echo "Processing chr${chr}..."

  if [ -n "${PLINK2_ENV}" ]; then
    conda activate "${PLINK2_ENV}"
  fi

  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
    "${PLINK2_BIN}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
    "${PLINK2_BIN}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  "${PLINK2_BIN}" \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/heartfail_filtered_chr${chr}_gwas_id_pwas.txt" \
    --keep "${EUR_SAMPLES}" \
    --make-bed \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_heartfail_pwas"

  "${PLINK2_BIN}" \
    --bfile "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_heartfail_pwas" \
    --export A \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_012_heartfail_pwas"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."

#!/bin/bash

# Build 1000 Genomes EUR reference genotype files for the HERMES HFpEF PWAS
# SNP lists. Step 4 needs:
#   1. .bim files for the SNP order/IDs.
#   2. .raw additive genotype dosage files from PLINK2 --export A.
# Run after 2_HERMES_HFpEF_prepare_data_pwas.R.

set -euo pipefail

PAPER_DIR="/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
DATA_DIR="${PAPER_DIR}/Data/plink_snplist_by_gene"
WORKSPACE_DIR="${PAPER_DIR}/Results/hermes_hfpef_pwas/hermes_hfpef_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
GWAS_ID_DIR="${WORKSPACE_DIR}/hermes_hfpef_filtered_id"
EUR_SAMPLES="${PAPER_DIR}/Data/1000Genome/eur_ids.txt"
CHR_LIST="$(echo {1..22})"
PLINK2_BIN="/opt/anaconda3/envs/plink2/bin/plink2"

for chr in ${CHR_LIST}; do
  echo "Processing chr${chr}..."

  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
    "${PLINK2_BIN}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
    "${PLINK2_BIN}" --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  # Keep only SNPs used by the matched weights/GWAS input and only EUR samples.
  "${PLINK2_BIN}" \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/hermes_hfpef_filtered_chr${chr}_gwas_id_pwas.txt" \
    --keep "${EUR_SAMPLES}" \
    --make-bed \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_hermes_hfpef_pwas"

  # Export additive 0/1/2 dosages. The filename keeps the historical "_012"
  # suffix because step 4 already reads this path.
  "${PLINK2_BIN}" \
    --bfile "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_hermes_hfpef_pwas" \
    --export A \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_012_hermes_hfpef_pwas"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."

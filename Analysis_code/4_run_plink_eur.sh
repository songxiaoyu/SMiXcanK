#!/bin/bash

# Explicitly source conda setup
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate plink2

# Base paths
CHR_LIST=$(echo {1..22})
DATA_DIR="/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene"
GWAS_ID_DIR="/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi_alpha08/bcac2020_filtered_id"

# Path to your EUR ID list (Make sure you created this file in Step 1)
# The file should contain one Sample ID per line

EUR_samples="/Users/zhusinan/Downloads/Data2/1000Genome/eur_ids.txt"

for chr in $CHR_LIST; do
  conda activate plink2
  echo "Processing chr${chr}..."

  # Decompress files if needed (only once per chr)
  # Added -f to force overwrite if file exists to prevent errors on re-runs
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pgen.zst" ]; then
      plink2 --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pgen.zst" "${DATA_DIR}/chr${chr}_hg38.pgen"
  fi
  if [ -f "${DATA_DIR}/chr${chr}_hg38.pvar.zst" ]; then
      plink2 --zst-decompress "${DATA_DIR}/chr${chr}_hg38.pvar.zst" "${DATA_DIR}/chr${chr}_hg38.pvar"
  fi

  # Extract SNPs AND Filter for EUR Samples
  # Added --keep "${EUR_samples}"
  plink2 \
    --pfile "${DATA_DIR}/chr${chr}_hg38" \
    --extract "${GWAS_ID_DIR}/bcac2020_filtered_chr${chr}_gwas_id_pi3_alpha08.txt" \
    --keep "${EUR_samples}" \
    --make-bed \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_pi3_alpha08_eur"

  # Convert to 012 format (additive coding)
  # Note: The input file here is already filtered for EUR, so no extra flag needed
  conda activate plink
  plink \
    --bfile "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_pi3_alpha08_eur" \
    --recodeA \
    --out "${GWAS_ID_DIR}/filtered_chr${chr}_hg38_012_pi3_alpha08_eur"

  echo "Finished chr${chr}."
done

echo "All chromosomes processed."
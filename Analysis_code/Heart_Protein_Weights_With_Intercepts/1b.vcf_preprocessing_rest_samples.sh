#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PLINK="${PLINK:-$HOME/Documents/plink2_mac_arm64_20260504/plink2}"
VCF="${VCF:-$HOME/Documents/vcf/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz}"
if [ ! -s "$VCF" ]; then
  VCF="$HOME/Documents/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
fi

SNPLIST="${SNPLIST:-$SCRIPT_DIR/sample_lists_heart_left_ventricle/heart_weights_pruned_snplist_100kb_1_r2_0.99.txt}"
SAMPLE_LIST="${SAMPLE_LIST:-$SCRIPT_DIR/sample_lists_heart_left_ventricle/heart_lv_vcf_sample_ids_not_in_original_EA110_and_in_vcf_351.txt}"
OUTDIR="${OUTDIR:-$SCRIPT_DIR/by_chr_heart_left_ventricle_351_samples_variant_id}"
BCFTOOLS_THREADS="${BCFTOOLS_THREADS:-4}"
PLINK_THREADS="${PLINK_THREADS:-8}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"
VARIANT_ID_TEMPLATE="${VARIANT_ID_TEMPLATE:-@:#:\$r:\$a}"

mkdir -p "$OUTDIR"

check_nonempty() {
  if [ ! -s "$1" ]; then
    echo "ERROR: missing or empty file: $1"
    exit 1
  fi
}

check_nonempty "$VCF"
check_nonempty "$SNPLIST"
check_nonempty "$SAMPLE_LIST"
if [ ! -x "$PLINK" ]; then
  echo "ERROR: PLINK2 executable not found or not executable: $PLINK"
  exit 1
fi

if [ ! -s "${VCF}.tbi" ] && [ ! -s "${VCF}.csi" ]; then
  echo "Indexing input VCF..."
  bcftools index -t "$VCF"
fi

# bcftools expression parser does not handle paths with spaces after ID=@ well.
# Use a temporary no-space symlink to the SNP list.
TMP_SNPLIST="${TMPDIR:-/tmp}/heart_weights_pruned_snplist_100kb_1_r2_0.99.$$.$RANDOM.txt"
ln -sf "$SNPLIST" "$TMP_SNPLIST"
trap 'rm -f "$TMP_SNPLIST"' EXIT

echo "Using 351 Heart Left Ventricle VCF samples not in original EA110:"
echo "  $SAMPLE_LIST"
echo "Sample count: $(wc -l < "$SAMPLE_LIST")"
echo "SNP list: $SNPLIST"
echo "SNP count: $(wc -l < "$SNPLIST")"
echo "Output directory: $OUTDIR"
echo "bcftools threads: $BCFTOOLS_THREADS"
echo "PLINK threads: $PLINK_THREADS"
echo "Keep intermediate VCF/PGEN files: $KEEP_INTERMEDIATE"
echo "PLINK variant ID template for dosage columns: $VARIANT_ID_TEMPLATE"
echo "Note: PLINK2 cannot read streamed VCF here, so a reduced temporary VCF is created per chromosome and deleted after dosage export unless KEEP_INTERMEDIATE=1."

for CHR in {1..22}; do
  echo "===== Processing Heart Left Ventricle 351 samples chr${CHR} ====="

  FILTERED_VCF="$OUTDIR/GTEx_HLV351_chr${CHR}_heartweights_snps_tmp.vcf.gz"
  FILTERED_TBI="$FILTERED_VCF.tbi"
  NOMISS_PREFIX="$OUTDIR/chr${CHR}_HLV351_nomiss"
  NOMISS_PGEN="${NOMISS_PREFIX}.pgen"
  DOSAGE_PREFIX="$OUTDIR/chr${CHR}_HLV351_dosage_nomiss"
  DOSAGE_RAW="${DOSAGE_PREFIX}.raw"

  # 1. Extract chromosome + 351 samples + heart-weights SNP list + SNP filters.
  # This temporary VCF is much smaller than the full chr/sample subset.
  if [ -s "$FILTERED_VCF" ]; then
    echo "Step 1 exists, skipping: $FILTERED_VCF"
  else
    echo "Step 1: creating reduced temporary VCF chr${CHR}"
    bcftools view \
      --threads "$BCFTOOLS_THREADS" \
      -r chr${CHR} \
      -S "$SAMPLE_LIST" \
      -I \
      -m2 -M2 -v snps \
      -i "ID=@$TMP_SNPLIST" \
      -Ou \
      "$VCF" \
    | bcftools view \
      --threads "$BCFTOOLS_THREADS" \
      -i '(STRLEN(REF)==1) && (STRLEN(ALT)==1) && (ALT!="A" || REF!="T") && (ALT!="T" || REF!="A") && (ALT!="C" || REF!="G") && (ALT!="G" || REF!="C")' \
      -Oz \
      -o "$FILTERED_VCF"
  fi
  check_nonempty "$FILTERED_VCF"

  # 2. Index temporary VCF.
  if [ -s "$FILTERED_TBI" ]; then
    echo "Step 2 exists, skipping: $FILTERED_TBI"
  else
    echo "Step 2: indexing reduced temporary VCF chr${CHR}"
    bcftools index --threads "$BCFTOOLS_THREADS" -t "$FILTERED_VCF"
  fi
  check_nonempty "$FILTERED_TBI"

  # 3. Remove SNPs with ANY missing genotype in this 351-sample set.
  if [ -s "$NOMISS_PGEN" ]; then
    echo "Step 3 exists, skipping: $NOMISS_PGEN"
  else
    echo "Step 3: removing missing SNPs chr${CHR}"
    "$PLINK" \
      --threads "$PLINK_THREADS" \
      --vcf "$FILTERED_VCF" \
      --geno 0 \
      --set-all-var-ids "$VARIANT_ID_TEMPLATE" \
      --make-pgen \
      --out "$NOMISS_PREFIX"
  fi
  check_nonempty "$NOMISS_PGEN"

  # 4. Export dosage/additive genotype matrix after missingness filtering.
  if [ -s "$DOSAGE_RAW" ]; then
    echo "Step 4 exists, skipping: $DOSAGE_RAW"
  else
    echo "Step 4: exporting dosage chr${CHR} with variant-ID columns"
    "$PLINK" \
      --threads "$PLINK_THREADS" \
      --pfile "$NOMISS_PREFIX" \
      --export A \
      --out "$DOSAGE_PREFIX"
  fi
  check_nonempty "$DOSAGE_RAW"

  if [ "$KEEP_INTERMEDIATE" = "0" ]; then
    echo "Cleaning intermediate files chr${CHR}"
    rm -f "$FILTERED_VCF" "$FILTERED_TBI" \
      "${NOMISS_PREFIX}.pgen" "${NOMISS_PREFIX}.pvar" "${NOMISS_PREFIX}.psam" "${NOMISS_PREFIX}.log"
  fi

  echo "===== Finished Heart Left Ventricle 351 samples chr${CHR} ====="
done

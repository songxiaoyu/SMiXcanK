# Heart Protein Weights With Intercepts

This folder contains the updated heart protein prediction workflow. It keeps the
old training setup, but saves intercepts and supports prediction in the 351 Heart
Left Ventricle samples.

## Main Inputs

- Original GTEx WGS VCF:
  `~/Documents/vcf/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz`
- 351 Heart Left Ventricle sample list:
  `New generated files/codes/sample_lists_heart_left_ventricle/heart_lv_vcf_sample_ids_not_in_original_EA110_and_in_vcf_351.txt`
- Heart-weight SNP list:
  `New generated files/codes/sample_lists_heart_left_ventricle/heart_weights_pruned_snplist_100kb_1_r2_0.99.txt`
- 351 sample cell fractions:
  `Heart/GTEx_Pi_Estimate/RNA_all_2celltypes_351_pi.tsv`

## Scripts

### `0_convert_HLV351_dosage_rsid_to_variant_id.R`

Converts already-created 351-sample dosage files from rsID columns to variant-ID
columns.

Use this only for old dosage files that look like:

```text
rs123_A
```

The converted files use the weight-compatible format:

```text
chr:pos:ref:alt_A
```

Current output folder:

```text
New generated files/codes/by_chr_heart_left_ventricle_351_samples_variant_id
```

### `1b.vcf_preprocessing_rest_samples.sh`

Generates dosage files for the 351 Heart Left Ventricle samples directly from the
VCF.

This script now sets PLINK variant IDs before exporting dosage:

```bash
--set-all-var-ids @:#:$r:$a
```

So future dosage files should already have columns like:

```text
22:17855114:C:T_C
```

This means the separate conversion script should not be needed after rerunning
`1b`.

### `2_train_heart_protein_weights_with_intercepts.R`

Retrains the heart protein models using the old workflow and saves intercepts.

It outputs:

- SMiXcan cell-type SNP weights
- intercepts for Tissue, Cardiomyocytes, and Other
- model terms
- diagnostics/skipped files
- PrediXcan-style tissue SNP weights from the Tissue model

It also includes the old combine logic. To combine per-chromosome outputs, run it
with:

```bash
PWAS_COMBINE_ONLY=1
```

### `run_train_heart_protein_weights_with_intercepts_by_chr.sh`

Runs script `2` chromosome-by-chromosome in parallel, then combines the
per-chromosome files using script `2` in combine-only mode.

Default training setup:

```text
Pruning: 100kb, r2 = 0.99
MiXcan alpha: 0.5
Lambda choice: lambda.min
Parallel jobs: 4
```

### `3_predict_celltype_heart_protein_HLV351_with_intercepts.R`

Predicts cell-type-specific protein levels in the 351 samples.

Formula:

```text
pred_cardiomyocytes = intercept_cardiomyocytes + sum(dosage * weight_cardiomyocytes)
pred_other          = intercept_other          + sum(dosage * weight_other)
```

Main output:

```text
Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts/predicted_heart_protein_HLV351_all_chr_with_intercepts.csv
```

### `run_predict_heart_protein_HLV351_with_intercepts.sh`

Small wrapper for script `3`.

### `4_construct_smixcan_tissue_level_prediction_HLV351.R`

Constructs tissue-level SMiXcan predictions by combining the cell-type-specific
predictions with the 351-sample PI estimates.

Formula:

```text
pred_tissue = pi_cardiomyocytes * pred_cardiomyocytes + pi_other * pred_other
```

Main output:

```text
Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts/predicted_heart_protein_HLV351_tissue_level_with_intercepts.csv
```

### `5_predict_predixcan_heart_protein_HLV351_with_intercepts.R`

Runs the PrediXcan-style tissue model directly in the 351 samples.

Formula:

```text
pred_tissue_predixcan = intercept_tissue + sum(dosage * weight_tissue)
```

Main output:

```text
Results/heart_protein_weights/predicted_protein_HLV351_predixcan_with_intercepts/predicted_heart_protein_HLV351_predixcan_all_chr_with_intercepts.csv
```

## Current Completed Outputs

### SMiXcan Cell-Type Prediction

```text
Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts/predicted_heart_protein_HLV351_all_chr_with_intercepts.csv
```

Rows: `2,122,848`

### SMiXcan Tissue-Level Prediction

```text
Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts/predicted_heart_protein_HLV351_tissue_level_with_intercepts.csv
```

Rows: `2,122,848`

### PrediXcan Tissue-Level Prediction

```text
Results/heart_protein_weights/predicted_protein_HLV351_predixcan_with_intercepts/predicted_heart_protein_HLV351_predixcan_all_chr_with_intercepts.csv
```

Rows: `2,122,848`

## Notes

- The old standalone combine script was merged into script `2` and archived in
  `Archive/`.
- The original old weights used `chr:pos:ref:alt`, not rsID.
- The 351 dosage files must use the same variant-ID format to match the weights.
- Intercepts are included in all prediction scripts.

# HERMES DCM PWAS Analysis

This folder contains scripts for running the PWAS/S-MiXcan workflow on the HERMES
DCM GWAS summary statistics.

## Input

Local HERMES folder:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Heart/HERMES
```

GWAS summary table:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR.tsv.gz
```

GWAS build:

```text
hg19 / GRCh37
```

Sample counts:

```text
cases    = 14,256
controls = 1,199,156
```

Relevant GWAS columns:

```text
chr, pos_b37, A1, A2, A1_beta, se, pval
```

`A1_beta` is interpreted as the effect size for `A1`. Downstream scripts use
`A1` as `Effect.Gwas` and `A2` as `Baseline.Gwas`.

## Outputs

Main output directory:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin
```

Figures:

```text
Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin
```

## Run The Workflow

Run from the repository root.

### Option A. Run one file at a time

#### 1. Liftover HERMES GWAS from hg19 to hg38

Input:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR.tsv.gz
Heart/Data/1000g_b38_snpIDs.txt
```

Script:

```bash
python3 Analysis_code/HERMES/1_liftover_hermes_gwas.py
```

Output:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR_hg38_rsid.tsv.gz
```

#### 2. Prepare HERMES PWAS input

Input:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR_hg38_rsid.tsv.gz
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Script:

```bash
Rscript Analysis_code/HERMES/2_HERMES_prepare_data_pwas.R
```

Output:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_input/chr{1..22}_mw_gwas_input_hermes_pwas.rds
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/hermes_filtered_chr{1..22}_gwas_id_pwas.txt
```

The RDS files contain the matched heart-protein weights and GWAS rows. The TXT
files contain SNP IDs used by step 3 to extract 1000 Genomes reference dosage.

#### 3. Extract 1000 Genomes EUR reference SNPs

Input:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/hermes_filtered_chr{1..22}_gwas_id_pwas.txt
Data/plink_snplist_by_gene/chr{1..22}_hg38.pgen/.pvar/.psam
Data/1000Genome/eur_ids.txt
```

Script:

```bash
bash Analysis_code/HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh
```

Output:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/filtered_chr{1..22}_hg38_hermes_pwas.bed/.bim/.fam
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/filtered_chr{1..22}_hg38_012_hermes_pwas.raw
```

The `.raw` files are additive 0/1/2 reference dosages used by step 4.

#### 4. Run HERMES PWAS/S-MiXcan association

Input:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_input/chr{1..22}_mw_gwas_input_hermes_pwas.rds
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/filtered_chr{1..22}_hg38_hermes_pwas.bim
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_filtered_id/filtered_chr{1..22}_hg38_012_hermes_pwas.raw
```

Script:

```bash
Rscript Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
```

Current settings:

```text
cases = 14256
controls = 1199156
GWAS family = binomial
regularization = fixed, scale 0.005
```

Output:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_chr{1..22}_result_pwas.csv
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas.csv
```

#### 5. Run Primo pattern annotation

Input:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas.csv
Data/ensembl38.txt
```

Script:

```bash
Rscript Analysis_code/HERMES/5_run_Primo_hermes_pwas.R
```

Output:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas_fixed_0p005_annotated.csv
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_table_pwas_fixed_0p005.csv
```

#### 6. Plot HERMES PWAS results

Input:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas.csv
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas_fixed_0p005_annotated.csv
```

Script:

```bash
Rscript Analysis_code/HERMES/6_plot_hermes_pwas.R
```

Output:

```text
Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_PWAS_fixed_0p005_QQ.pdf
Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_PWAS_fixed_0p005_QQ_split_venn.pdf
```

### Option B. Run the full workflow with one script

```bash
bash Analysis_code/HERMES/run_HERMES_analysis.sh
```

The current wrapper runs steps 4-6 only. Use it after steps 1-3 have already
created the lifted GWAS, matched input RDS files, and 1000 Genomes reference
dosage files.

## Notes

- This workflow uses these heart-protein weights:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

- To run one chromosome for testing:

```bash
Rscript Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R 1
```

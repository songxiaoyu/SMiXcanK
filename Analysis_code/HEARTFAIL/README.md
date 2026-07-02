# HEARTFAIL PWAS Analysis

This folder contains the HEARTFAIL GWAS version of the PWAS/S-MiXcan workflow.

## Input

Default GWAS file:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes.tsv.bgz
```

The raw GWAS columns are expected to include:

```text
variant, minor_allele, beta, se, pval
```

`variant` is parsed as `chr:pos:ref:alt` in hg19/GRCh37. The script uses
`minor_allele` as the GWAS effect allele, so if `minor_allele == REF`, the
baseline allele is `ALT`; if `minor_allele == ALT`, the baseline allele is
`REF`.

## Output

Default workspace:

```text
Results/heartfail_pwas/heartfail_workspace
```

Main outputs:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_table_pwas.csv
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ.pdf
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ_split_venn.pdf
```

## Run The Workflow

Run from the repository root.

### Option A. Run one file at a time

#### 1. Lift HEARTFAIL GWAS from hg19 to hg38

Input:

```text
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes.tsv.bgz
Heart/Data/1000g_b38_snpIDs.txt
```

Script:

```bash
python3 Analysis_code/HEARTFAIL/1_liftover_heartfail_gwas.py
```

Output:

```text
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes_hg38_rsid_chr_pos.txt
```

#### 2. Merge PWAS weights with HEARTFAIL GWAS

Input:

```text
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Script:

```bash
Rscript Analysis_code/HEARTFAIL/2_HEARTFAIL_prepare_data_pwas.R
```

Output:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_input/chr{1..22}_mw_gwas_input_heartfail_pwas.rds
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/heartfail_filtered_chr{1..22}_gwas_id_pwas.txt
```

#### 3. Build EUR LD reference files

Input:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/heartfail_filtered_chr{1..22}_gwas_id_pwas.txt
Data/plink_snplist_by_gene/chr{1..22}_hg38.pgen/.pvar/.psam
Data/1000Genome/eur_ids.txt
```

Script:

```bash
PLINK2_BIN=/opt/anaconda3/envs/plink2/bin/plink2 \
bash Analysis_code/HEARTFAIL/3_1000Genome_keep_eur_plink_heartfail_pwas.sh
```

Output:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/filtered_chr{1..22}_hg38_heartfail_pwas.bed/.bim/.fam
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/filtered_chr{1..22}_hg38_012_heartfail_pwas.raw
```

#### 4. Run S-MiXcan association

Input:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_input/chr{1..22}_mw_gwas_input_heartfail_pwas.rds
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/filtered_chr{1..22}_hg38_heartfail_pwas.bim
Results/heartfail_pwas/heartfail_workspace/heartfail_filtered_id/filtered_chr{1..22}_hg38_012_heartfail_pwas.raw
```

Script:

Use the current HEARTFAIL counts:

```bash
HEARTFAIL_N_CASES=1405 \
HEARTFAIL_N_CONTROLS=359789 \
Rscript Analysis_code/HEARTFAIL/4_HEARTFAIL_run_analysis_pwas.R
```

Output:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_chr{1..22}_result_pwas.csv
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
```

#### 5. Run Primo pattern inference

Input:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
Data/ensembl38.txt
```

Script:

```bash
Rscript Analysis_code/HEARTFAIL/5_run_Primo_heartfail_pwas.R
```

Output:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas_annotated.csv
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_table_pwas.csv
```

#### 6. Plot QQ and split Venn

Input:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas_annotated.csv
```

Script:

```bash
Rscript Analysis_code/HEARTFAIL/6_plot_heartfail_pwas.R
```

Output:

```text
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ.pdf
Figure/heartfail_pwas/HEARTFAIL_PWAS_QQ_split_venn.pdf
```

### Option B. Run the full workflow with one script

```bash
bash Analysis_code/HEARTFAIL/run_HEARTFAIL_analysis.sh
```

The wrapper uses:

```text
cases = 1405
controls = 359789
regularization = fixed, scale 0.1
parallel step 4 jobs = 4
```

## Notes

Confirm the HEARTFAIL case/control counts before rerunning if the GWAS input file changes.

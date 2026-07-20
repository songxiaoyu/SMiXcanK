# Analysis Pipeline Guide

This folder contains the analysis scripts used for the `SMiXcan` paper and
follow-up PWAS workflows. The code is organized into five active workflow
folders:

- `Breast_BCAC_2CellTypes/`: breast cancer risk / BCAC two-cell-type workflow.
- `Breast_BCAC_3CellTypes/`: breast cancer risk / BCAC three-cell-type workflow.
- `Heart_Protein_Weights/`: heart protein weight training workflow.
- `HERMES/`: HERMES DCM GWAS PWAS workflow.
- `HERMES_HF/`: HERMES HF GWAS PWAS workflow.

Older exploratory heart PWAS folders and superseded result/figure snapshots have
been moved outside the GitHub code folder to `Paper_SMiXcan/Archive`.

The scripts were originally written for a project directory with several large
reference datasets stored outside GitHub. This README explains the intended run
order, what each script does, and which files need to exist before each stage
can run.

## Before You Start

This repository is the code layer, not a full self-contained data bundle. The
active scripts expect large data and intermediate files under the external
`Paper_SMiXcan` project directory, usually controlled by:

```bash
PAPER_SMIXCAN_DIR=/path/to/Paper_SMiXcan
```

On the current analysis machine, the checked workflow directories are under the
Dropbox-backed `Paper_SMiXcan` folder.

## Expected Folder Layout

The current scripts assume the following project layout:

```text
Paper_SMiXcan/
├── Data/
│   ├── BCAC/
│   │   ├── chr1_icogs_hg38.csv
│   │   ├── ...
│   │   └── chr22_icogs_hg38.csv
│   ├── DRIVE/
│   ├── GTEx/
│   └── ensembl38.txt
├── Heart/
│   ├── Data/
│   ├── GTEx_Pi_Estimate/
│   └── HERMES/
├── Results/
│   ├── BayesDeBulk_pi_3ct_GTEx.tsv
│   ├── SMiXcanK_results/
│   ├── drive_result_full_lam_new.csv
│   ├── heart_protein_weights/
│   ├── hermes_pwas/
│   └── hermes_hf_pwas/
├── Figure/
└── Github/
    └── Analysis_code/
```

The active workflows currently use these external data/workspace directories:

- `Data/GTEx`
- `Data/BCAC`
- `Data/DRIVE`
- `Data/1000Genome`
- `Results/2pi_workspace`
- `Results/3pi_workspace`
- `Heart/GTEx_Pi_Estimate`
- `Heart/HERMES`
- `Heart/HERMES_HF`
- `Heart/Data`
- `New generated files/codes`
- `Results/heart_protein_weights`
- `Results/hermes_pwas`
- `Results/hermes_hf_pwas`

If you run the scripts on a different machine, set `PAPER_SMIXCAN_DIR` to the
local path of the external `Paper_SMiXcan` project directory before running the
workflow scripts.

## Package Setup

Install the package from the repository root before running the analysis scripts:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "ggplot2",
  "ggrepel", "cowplot", "ggforce", "bacon", "glmnet", "lme4",
  "doParallel", "doRNG", "janitor", "tibble", "DBI", "RSQLite",
  "MASS", "ACAT", "Primo"
))
install_github("petraf01/BayesDeBulk/BayesDeBulk")

install.packages(".", repos = NULL, type = "source")
```

```
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("limma")

install.packages("MASS")
install.packages("matrixStats")
install.packages("nnls")
install.packages("R.methodsS3")
install.packages("lcmix",repos="http://r-forge.r-project.org")

devtools::install_github("kjgleason/Primo")
library("Primo")
```

Many scripts use `library(SMiXcan)`, so the package should be installed before running them.

The active association function is the package-level `SMiXcan_assoc_test_K()`.
Do not source local copies of `R/S-MiXcan_K.R` from workflow scripts.

## Breast Cancer Risk / BCAC Two-Cell-Type Pipeline

The breast cancer risk two-cell analysis lives in
[Breast_BCAC_2CellTypes](./Breast_BCAC_2CellTypes).

### Overview

This workflow:

1. trains a two-cell-type prediction model in GTEx breast tissue
2. evaluates agreement between MiXcan and S-MiXcan in the DRIVE dataset
3. prepares BCAC summary statistics for chromosome-wise analysis
4. builds European-only LD reference subsets from 1000 Genomes
5. runs the S-MiXcan association scan across all chromosomes
6. annotates results and performs Primo pattern inference
7. generates the main figures

### Script Order

1. [Breast_BCAC_2CellTypes/1_pi_estimate_2.R](./Breast_BCAC_2CellTypes/1_pi_estimate_2.R)

This file is only a note. The two-cell proportions used here are derived in the training workflow from the 3-cell BayesDeBulk estimates, so there is no standalone 2-cell `pi` estimation script to run.

2. [Breast_BCAC_2CellTypes/2_train_model_2.R](./Breast_BCAC_2CellTypes/2_train_model_2.R)

Purpose:
- loads GTEx breast expression, covariates, genotype, and PredictDB elastic-net SNP sets
- loads `Results/BayesDeBulk_pi_3ct_GTEx.tsv`
- collapses the 3-cell proportions into epithelial versus other
- trains the two-cell MiXcan model gene by gene

Main output:
- `Paper_SMiXcan/Results/weights_miXcan_full_pi2.csv`

3. [Breast_BCAC_2CellTypes/3_DRIVE_Analysis_gtex.R](./Breast_BCAC_2CellTypes/3_DRIVE_Analysis_gtex.R)

Purpose:
- loads the trained two-cell weights
- loads DRIVE genotype and phenotype data
- runs single-SNP GWAS in DRIVE
- compares individual-level MiXcan and summary-statistics S-MiXcan

Main output:
- `Paper_SMiXcan/Results/drive_result_full_lam_new.csv`

Notes:
- this script contains a small local helper called `smixcan_assoc_for_drive()`
- that helper is only for this comparison script and is not part of the package API

4. [Breast_BCAC_2CellTypes/4_BCAC_prepare_data.R](./Breast_BCAC_2CellTypes/4_BCAC_prepare_data.R)

Purpose:
- reads chromosome-wise BCAC GWAS summary statistics
- matches BCAC SNPs to the trained two-cell weight table
- writes chromosome-specific merged inputs

Input location:
- `Paper_SMiXcan/Data/BCAC/chr*_icogs_hg38.csv`

Source:
- BCAC GWAS summary statistics were downloaded from the BCAC summary results page: `https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations`

Main outputs:
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_filtered_id/`
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_input/`

5. [Breast_BCAC_2CellTypes/5_1000Genome_keep_eur_plink.sh](./Breast_BCAC_2CellTypes/5_1000Genome_keep_eur_plink.sh)

Purpose:
- uses the SNP lists from step 4
- extracts European-only reference genotypes from 1000 Genomes with PLINK

Main output location:
- `Paper_SMiXcan/Results/2pi_workspace/bcac2020_filtered_id/`

This step requires PLINK and the 1000 Genomes reference files under
`PAPER_SMIXCAN_DIR`.

6. [Breast_BCAC_2CellTypes/6_BCAC2020_run_analysis_2.R](./Breast_BCAC_2CellTypes/6_BCAC2020_run_analysis_2.R)

Purpose:
- loads the chromosome-wise merged BCAC inputs
- loads the European-only LD reference panels from step 5
- runs `SMiXcan_assoc_test_K()` gene by gene across chromosomes
- combines all chromosome outputs

Main output location:
- `Paper_SMiXcan/Results/bcac2020_result/`

7. [Breast_BCAC_2CellTypes/7_run_Primo.R](./Breast_BCAC_2CellTypes/7_run_Primo.R)

Purpose:
- reads `bcac2020_result_pi2.csv`
- adds gene name and cytoband annotation using `Data/ensembl38.txt`
- runs `infer_celltype_patterns()`
- prepares the supplementary table for the two-cell analysis

Main outputs in `Paper_SMiXcan/Results/SMiXcanK_results`:
- `bcac2020_result_pi2_annotated.csv`
- `tableS2.csv`

8. [Breast_BCAC_2CellTypes/8_plot.R](./Breast_BCAC_2CellTypes/8_plot.R)

Purpose:
- uses the annotated two-cell BCAC output and the DRIVE comparison output
- generates the main manuscript figures

Main outputs in `Paper_SMiXcan/Figure`:
- `Figure1_ABC_scatter.pdf`
- `Figure2_BC.pdf`
- `Final_Figure_ABC.pdf`

## Breast Cancer Risk / BCAC Three-Cell-Type Pipeline

The breast cancer risk three-cell analysis lives in
[Breast_BCAC_3CellTypes](./Breast_BCAC_3CellTypes).

### Overview

This workflow:

1. estimates three-cell-type proportions in GTEx breast tissue
2. trains a three-cell MiXcan model
3. validates three-cell MiXcan and S-MiXcan agreement in DRIVE
4. prepares BCAC chromosome-wise inputs
5. builds European-only LD reference subsets
6. runs the three-cell S-MiXcan scan
7. annotates the results and runs Primo
8. generates the supplementary figures

### Script Order

1. [Breast_BCAC_3CellTypes/1_pi_estimation_gtex_3.R](./Breast_BCAC_3CellTypes/1_pi_estimation_gtex_3.R)

Purpose:
- estimates cell fractions with `pi_estimation_K()`
- uses a marker list for adipose/endothelial, fibroblast, and epithelial compartments

Main output:
- `Paper_SMiXcan/Results/BayesDeBulk_pi_3ct_GTEx.tsv`

2. [Breast_BCAC_3CellTypes/2_train_model_3.R](./Breast_BCAC_3CellTypes/2_train_model_3.R)

Purpose:
- loads GTEx expression, covariates, genotype, and PredictDB elastic-net SNP sets
- loads the three-cell BayesDeBulk proportions
- trains a three-cell MiXcan model gene by gene

Main output:
- `Paper_SMiXcan/Results/weights_miXcan_full_pi3.csv`

DRIVE validation: [Breast_BCAC_3CellTypes/3_DRIVE_Analysis_gtex_3.R](./Breast_BCAC_3CellTypes/3_DRIVE_Analysis_gtex_3.R)

Purpose:
- runs a single-SNP GWAS in the individual-level DRIVE chromosome 21 data
- compares three-cell individual-level MiXcan with S-MiXcan
- allele-aligns DRIVE dosages to the three-cell prediction weights

Main output:
- `Paper_SMiXcan/Results/drive_result_pi3.csv`

3. [Breast_BCAC_3CellTypes/3_BCAC2020_prepare_data_3.R](./Breast_BCAC_3CellTypes/3_BCAC2020_prepare_data_3.R)

Purpose:
- reads chromosome-wise BCAC GWAS summary statistics
- matches BCAC SNPs to the three-cell weight table
- writes chromosome-specific merged inputs

Input location:
- `Paper_SMiXcan/Data/BCAC/chr*_icogs_hg38.csv`

Source:
- BCAC GWAS summary statistics were downloaded from the BCAC summary results page: `https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-associations`

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/`

Note:
- these archived Dropbox workspace folders contain the available 3-cell chromosome-specific merged inputs and SNP ID files
- the original workspace folder keeps the name `baca2020_input/` in Dropbox because that is how it was named in the local run directory

4. [Breast_BCAC_3CellTypes/4_run_plink_eur_3.sh](./Breast_BCAC_3CellTypes/4_run_plink_eur_3.sh)

Purpose:
- uses the SNP lists from step 3
- extracts the 1000 Genomes reference genotypes used by the kept `3pi` workflow with PLINK

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/bcac2020_filtered_id/`

5. [Breast_BCAC_3CellTypes/5_BCAC2020_run_analysis_3.R](./Breast_BCAC_3CellTypes/5_BCAC2020_run_analysis_3.R)

Purpose:
- loads the chromosome-wise merged BCAC inputs
- loads the European-only LD reference panels from step 4
- runs `SMiXcan_assoc_test_K()` gene by gene across chromosomes
- combines all chromosome outputs

Main output locations:
- `Paper_SMiXcan/Results/3pi_workspace/bcac2020_result/`

6. [Breast_BCAC_3CellTypes/6_BCACresult.R](./Breast_BCAC_3CellTypes/6_BCACresult.R)

Purpose:
- reads the combined three-cell result table
- adds gene and cytoband annotation
- runs `infer_celltype_patterns()`
- prepares the supplementary table for the three-cell analysis

Main outputs in `Paper_SMiXcan/Results/SMiXcanK_results`:
- `bcac2020_result_pi3_annotated.csv`
- `tableS3.csv`

7. [Breast_BCAC_3CellTypes/7_plotS1.R](./Breast_BCAC_3CellTypes/7_plotS1.R)

Purpose:
- uses the annotated three-cell BCAC output and the DRIVE comparison output
- generates the supplementary figures

Main outputs in `Paper_SMiXcan/Figure`:
- `FigureS1_111.pdf`
- `FigureS2_111.pdf`

## Heart Protein PWAS Workflows

The heart protein workflows use the current two-component setup:

```text
cardiomyocytes
other = fibroblasts + endothelial + smooth muscle/pericyte + immune
```

The active workflow folders are:

```text
Heart_Protein_Weights/
HERMES/
HERMES_HF/
```

### Heart Protein Weight Training

See [Heart_Protein_Weights/README.md](./Heart_Protein_Weights/README.md).

Main scripts:

1. [Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh](./Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh)
2. [Heart_Protein_Weights/2_train_heart_protein_weights.R](./Heart_Protein_Weights/2_train_heart_protein_weights.R)
3. [Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R](./Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R)

The current active weight source is documented in:

```text
Heart_Protein_Weights/CURRENT_WEIGHTS_SOURCE.md
```

### HERMES DCM PWAS

See [HERMES/README.md](./HERMES/README.md).

Main scripts:

1. [HERMES/1_liftover_hermes_gwas.py](./HERMES/1_liftover_hermes_gwas.py)
2. [HERMES/2_HERMES_prepare_data_pwas.R](./HERMES/2_HERMES_prepare_data_pwas.R)
3. [HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh](./HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh)
4. [HERMES/4_HERMES_run_analysis_pwas.R](./HERMES/4_HERMES_run_analysis_pwas.R)
5. [HERMES/5_run_Primo_hermes_pwas.R](./HERMES/5_run_Primo_hermes_pwas.R)
6. [HERMES/6_plot_hermes_pwas.R](./HERMES/6_plot_hermes_pwas.R)

Convenience runner:

```bash
bash Analysis_code/HERMES/run_HERMES_analysis.sh
```

HERMES default sample counts:

```text
cases    = 14,256
controls = 1,199,156
```

Step 4 supports:

```bash
HERMES_ASSOC_REG_MODE=fixed
HERMES_ASSOC_REG_MODE=estimate
```

### HERMES HF PWAS

See [HERMES_HF/README.md](./HERMES_HF/README.md).

Main scripts:

1. [HERMES_HF/1_liftover_hermes_hf_gwas.py](./HERMES_HF/1_liftover_hermes_hf_gwas.py)
2. [HERMES_HF/2_HERMES_HF_prepare_data_pwas.R](./HERMES_HF/2_HERMES_HF_prepare_data_pwas.R)
3. [HERMES_HF/3_1000Genome_keep_eur_plink_hermes_hf_pwas.sh](./HERMES_HF/3_1000Genome_keep_eur_plink_hermes_hf_pwas.sh)
4. [HERMES_HF/4_HERMES_HF_run_analysis_pwas.R](./HERMES_HF/4_HERMES_HF_run_analysis_pwas.R)
5. [HERMES_HF/5_run_Primo_hermes_hf_pwas.R](./HERMES_HF/5_run_Primo_hermes_hf_pwas.R)
6. [HERMES_HF/6_plot_hermes_hf_pwas.R](./HERMES_HF/6_plot_hermes_hf_pwas.R)

HERMES HF default sample counts:

```text
cases    = 132,176
controls = 1,553,537
```

Step 4 supports:

```bash
HERMES_HF_ASSOC_REG_MODE=fixed
HERMES_HF_ASSOC_REG_MODE=estimate
```

## Starting a New Heart PWAS Analysis From Scratch

This section is for a user who wants to start from raw/new inputs, train the
heart protein weights, and then run PWAS on a GWAS summary table.

### Step 0. Decide the analysis root

All heart PWAS scripts assume a project directory named by `PAPER_SMIXCAN_DIR`.
Set it before running anything:

```bash
export PAPER_SMIXCAN_DIR=/path/to/Paper_SMiXcan
```

The expected data layout is:

```text
${PAPER_SMIXCAN_DIR}/
  Data/
    1000Genome/
      eur_ids.txt
    ensembl38.txt
  Heart/
    GTEx_Pi_Estimate/
      Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData
      BayesDeBulk_pi.tsv
      GTEx_heart_EA_subject_ids.txt
    Data/
      <GWAS files>
    HERMES/
      <HERMES files, if applicable>
  New generated files/
    covariate_EA_with_age.txt
    codes/
      by_chr/
      by_chr_nomiss/
      pruned_by_chr/
  Results/
  Figure/
```

### Step 1. Prepare genotype dosage for training

The training script expects per-chromosome GTEx dosage `.raw` files and matching
variant annotation files. The current moderate-pruning workflow is:

```bash
bash Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh
```

Things to check or modify:

- `PAPER_SMIXCAN_DIR`
- `PWAS_PLINK2_ENV`, `PLINK2_BIN`, or `PLINK_BIN` if PLINK is not on `PATH`
- genotype input location under `New generated files/codes`
- pruning settings inside the script if you want a different SNP density

Current preferred setting:

```text
100kb window
r2 = 0.99
alpha = 0.5
lambda = lambda.min
```

### Step 2. Train heart protein weights

Main script:

```bash
bash Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh
```

or run a single chromosome:

```bash
PWAS_CHR_FILTER=1 \
PWAS_OUTPUT_SUFFIX=_moderate_100kb_r2_0.99_alpha0.5_lambdamin \
PWAS_LAMBDA_CHOICE=lambda.min \
PWAS_MIXCAN_ALPHA=0.5 \
Rscript Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R
```

Important inputs in `2_train_heart_protein_weights.R`:

```text
protein_file:
  Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData

pi_file:
  Heart/GTEx_Pi_Estimate/BayesDeBulk_pi.tsv

covariate_file:
  New generated files/covariate_EA_with_age.txt

ea_keep_file:
  Heart/GTEx_Pi_Estimate/GTEx_heart_EA_subject_ids.txt

ensembl_file:
  Data/ensembl38.txt

rsid_annotation_file:
  Data/plink_snplist_by_gene/ALL_genes_snplist_with_rsids.csv

geno_raw_dir:
  New generated files/codes/by_chr_nomiss
```

Things a new user may need to change:

- input protein RData filename
- covariate filename and covariate columns
- pi/cell-type proportion file
- dosage folder
- rsID annotation file
- `PWAS_CHR_FILTER` for test runs
- `PWAS_OUTPUT_SUFFIX` so new experiments do not overwrite old results
- `PWAS_LAMBDA_CHOICE` and `PWAS_MIXCAN_ALPHA`

Current output directory:

```text
Results/heart_protein_weights/training_model_weights/
```

Per-chromosome outputs are named like:

```text
weights_heart_protein_cardiomyocytes_other_chr<CHR><SUFFIX>.csv
skipped_heart_protein_cardiomyocytes_other_chr<CHR><SUFFIX>.csv
diagnostics_heart_protein_cardiomyocytes_other_chr<CHR><SUFFIX>.csv
```

### Step 3. Combine chromosome-level weights

After all chromosomes finish:

```bash
Rscript Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R
```

Current combined weight file:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

If you use a new suffix, update `PWAS_OUTPUT_PREFIX` or the combine script so it
points to the correct per-chromosome files.

### Step 4. Prepare a GWAS for PWAS

For a new GWAS, create or copy a workflow folder under `Analysis_code/`, using
`HERMES/` or `HERMES_HF/` as a template.

You usually need these scripts:

```text
1_liftover_<trait>_gwas.py
2_<TRAIT>_prepare_data_pwas.R
3_1000Genome_keep_eur_plink_<trait>_pwas.sh
4_<TRAIT>_run_analysis_pwas.R
5_run_Primo_<trait>_pwas.R
6_plot_<trait>_pwas.R
```

Files/information to update:

```text
GWAS input file
GWAS genome build: hg19 or hg38
GWAS column names: chr, pos, effect allele, non-effect allele, beta, se, p
case count
control count
output workspace name
plot output folder
weight file path
```

If the GWAS is hg19, run liftover first. If it is already hg38 with rsIDs, step
1 can be skipped or simplified.

Minimum edit checklist for a new GWAS workflow:

```text
1_liftover_<trait>_gwas.py
  input GWAS file
  output lifted/annotated GWAS file
  source genome build and target genome build

2_<TRAIT>_prepare_data_pwas.R
  gwas_file
  weights_file
  workspace_dir
  GWAS column mapping
  allele harmonization columns

3_1000Genome_keep_eur_plink_<trait>_pwas.sh
  WORKSPACE_DIR
  DATA_DIR / 1000 Genomes PLINK reference path
  EUR sample list
  PLINK binary path or conda environment

4_<TRAIT>_run_analysis_pwas.R
  workspace_dir
  result_dir
  case count
  control count
  association regularization mode

5_run_Primo_<trait>_pwas.R
  result file
  annotation file
  PRIMO output file

6_plot_<trait>_pwas.R
  result file
  plot output directory
  plot title / trait label
```

### Step 5. Merge weights with GWAS

Example for HERMES:

```bash
Rscript Analysis_code/HERMES/2_HERMES_prepare_data_pwas.R
```

Example for HERMES HF:

```bash
Rscript Analysis_code/HERMES_HF/2_HERMES_HF_prepare_data_pwas.R
```

Things to check in the prepare script:

- `gwas_file`
- `weights_file`
- `workspace_dir`
- allele columns used for harmonization
- chromosome list

The prepare step writes:

```text
<workspace>/<trait>_input/chr*_mw_gwas_input_<trait>_pwas.rds
<workspace>/<trait>_filtered_id/*_gwas_id_pwas.txt
```

### Step 6. Build LD reference inputs

Example:

```bash
bash Analysis_code/HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh
```

or:

```bash
bash Analysis_code/HERMES_HF/3_1000Genome_keep_eur_plink_hermes_hf_pwas.sh
```

Things to check:

- `PAPER_SMIXCAN_DIR`
- `DATA_DIR`, usually `Data/plink_snplist_by_gene`
- `EUR_SAMPLES`, usually `Data/1000Genome/eur_ids.txt`
- `PLINK2_BIN` / `PLINK_BIN`
- `HERMES_WORKSPACE_DIR` or `HERMES_HF_WORKSPACE_DIR`

This step writes `.bim` and `.raw` files used by step 4.

### Step 7. Run association

Example HERMES:

```bash
HERMES_ASSOC_REG_MODE=estimate \
Rscript Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
```

Example HERMES HF:

```bash
HERMES_HF_ASSOC_REG_MODE=estimate \
Rscript Analysis_code/HERMES_HF/4_HERMES_HF_run_analysis_pwas.R
```

Association options:

```text
*_ASSOC_REG_MODE=fixed
*_ASSOC_REG_SCALE=0.1
```

or:

```text
*_ASSOC_REG_MODE=estimate
```

`estimate` uses the package defaults inside `SMiXcan_assoc_test_K()`.

For long runs, use the parallel runner:

```bash
bash Analysis_code/HERMES/run_HERMES_step4_parallel.sh
bash Analysis_code/HERMES_HF/run_HERMES_HF_step4_parallel.sh
```

### Step 8. Run PRIMO and plot

Example:

```bash
Rscript Analysis_code/HERMES/5_run_Primo_hermes_pwas.R
Rscript Analysis_code/HERMES/6_plot_hermes_pwas.R
```

or:

```bash
Rscript Analysis_code/HERMES_HF/5_run_Primo_hermes_hf_pwas.R
Rscript Analysis_code/HERMES_HF/6_plot_hermes_hf_pwas.R
```

Final outputs are written under:

```text
Results/<trait>_pwas/<workspace>/<trait>_result/
Figure/<trait>_pwas/
```

### Step 9. Quality checks before interpreting results

Before interpreting significant genes, check:

- number of tested genes/proteins
- `input_snp_num` distribution
- minimum p-value
- FDR < 0.10 count
- QQ plot
- whether top genes are driven by 1 SNP
- whether top genes have reasonable CV R2 in the training diagnostics
- whether GWAS allele harmonization and genome build are correct

For heart PWAS analyses, report FDR < 0.10 as the main significance threshold
unless the project specifies otherwise.

## Typical Run Order for Reproducing the Main Results

If you want the original breast cancer risk / BCAC paper outputs rather than every intermediate step, the practical run order is:

1. install the `SMiXcan` package from this repository
2. make sure Dropbox `Data`, `Results`, and `Figure` folders exist
3. make sure the external GTEx, BCAC, and 1000 Genomes resources are available under `PAPER_SMIXCAN_DIR`
4. run `Breast_BCAC_2CellTypes/2_train_model_2.R`
5. run `Breast_BCAC_2CellTypes/3_DRIVE_Analysis_gtex.R`
6. run `Breast_BCAC_2CellTypes/4_BCAC_prepare_data.R`
7. run `Breast_BCAC_2CellTypes/5_1000Genome_keep_eur_plink.sh`
8. run `Breast_BCAC_2CellTypes/6_BCAC2020_run_analysis_2.R`
9. run `Breast_BCAC_2CellTypes/7_run_Primo.R`
10. run `Breast_BCAC_2CellTypes/8_plot.R`

For the three-cell supplementary analysis, run the corresponding `Breast_BCAC_3CellTypes` scripts in numeric order.

For the heart PWAS analyses, train or confirm the heart protein weights first,
then run either the `HERMES` or `HERMES_HF` folder in numeric order.

## Important Path Caveats

There are still a few path assumptions to be aware of:

- Most scripts default to `PAPER_SMIXCAN_DIR=/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan`, but this can be overridden with the `PAPER_SMIXCAN_DIR` environment variable.
- The two-cell pipeline uses `Results/2pi_workspace` for chromosome-specific intermediate files.
- The three-cell pipeline uses `Results/3pi_workspace` for chromosome-specific intermediate files.
- The archived 3-cell workspace retains one historical folder typo, `baca2020_input/`, and the scripts follow that name so they match the stored run artifacts.

If your machine uses a different external project directory, set
`PAPER_SMIXCAN_DIR` consistently before running.

## Notes for New Users

- The scripts are intended to be run one file at a time, not as a single automated pipeline
- Most intermediate files are chromosome-specific and are written into the `Results/2pi_workspace`, `Results/3pi_workspace`, HERMES, or HERMES HF workspace directories
- Final cleaned outputs are written into Dropbox `Results`
- Final figures are written into Dropbox `Figure`
- The file [Document for downloading reference genome.docx](./Document%20for%20downloading%20reference%20genome.docx) contains the download/setup notes for the external reference resources used in this project, including the GTEx-related local data setup and the reference genome resources used by the PLINK steps
- The local research report `PWAS_RESEARCH_REPORT.md` is a diagnostic handoff document and is not part of the public package build

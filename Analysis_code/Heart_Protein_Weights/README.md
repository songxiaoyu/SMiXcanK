# Heart Protein Weights Training

This folder trains GTEx heart protein prediction weights for downstream
PWAS/S-MiXcan analyses, including HERMES DCM, HERMES HFpEF, and HEARTFAIL.

The retained model predicts heart protein abundance using two components:

1. Cardiomyocytes
2. Other cell types combined

The combined output weights are saved under:

```text
Results/heart_protein_weights/training_model_weights
```

## Current Retained Weights

Current retained weights file:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Current skipped-protein log:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_skipped_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Older per-chromosome, diagnostic, and superseded weight outputs were moved to:

```text
Archive/experiment/heart_protein_weights_training_model_weights_archive_20260702
```

Training settings:

```text
cell components = cardiomyocytes + other
genotype input  = GTEx EA heart dosage after moderate LD pruning
pruning         = 100 kb window, step 1, r2 = 0.99
MiXcan alpha    = 0.5
glmnet lambda   = lambda.min
```

See [CURRENT_WEIGHTS_SOURCE.md](./CURRENT_WEIGHTS_SOURCE.md) for the retained
weight source record.

## Inputs

Protein phenotype:

```text
Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData
```

Cell fraction estimates:

```text
Heart/GTEx_Pi_Estimate/BayesDeBulk_pi.tsv
```

Covariates:

```text
New generated files/covariate_EA_with_age.txt
```

EA sample list:

```text
Heart/GTEx_Pi_Estimate/GTEx_heart_EA_subject_ids.txt
```

Gene annotation:

```text
Data/ensembl38.txt
```

No-missing GTEx dosage source:

```text
New generated files/codes/by_chr_nomiss
```

Moderate-pruned dosage used by the retained model:

```text
New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99
```

## Scripts

1. `1_prune_moderate_heart_protein_dosage.sh`
   - Creates moderately LD-pruned dosage files from no-missing GTEx dosage.
   - The script is generic; pass pruning parameters when reproducing the
     retained model.

2. `2_train_heart_protein_weights.R`
   - Trains chromosome-level MiXcan protein prediction weights.
   - Uses the two cell components, covariates, protein phenotype, and genotype
     dosage.

3. `run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh`
   - Runs `2_train_heart_protein_weights.R` across chromosomes with the retained
     model settings.

4. `3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R`
   - Combines per-chromosome weights into one full weights table.

5. `run_heart_protein_weights_pipeline.sh`
   - Runs pruning, chromosome-level training, and combining in one command.

## Run The Workflow

Run from the repository root.

### Option A. Run one file at a time

Step 1: create the retained moderate-pruned dosage files:

```bash
PWAS_PRUNE_KB=100 \
PWAS_PRUNE_STEP=1 \
PWAS_PRUNE_R2=0.99 \
bash Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh
```

Step 2: train chromosome-level weights:

```bash
bash Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh
```

Step 3: combine chromosome-level weights:

```bash
Rscript Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R
```

### Option B. Run the full workflow with one script

```bash
bash Analysis_code/Heart_Protein_Weights/run_heart_protein_weights_pipeline.sh
```

The wrapper uses the retained pruning settings:

```text
PWAS_PRUNE_KB=100
PWAS_PRUNE_STEP=1
PWAS_PRUNE_R2=0.99
```

## Output Format

The final weights table contains SNP-level nonzero model weights. Important
columns include:

```text
gene_id
gene_name
varID
chr
pos
ref_allele
eff_allele
weight_cardiomyocytes
weight_other
type
```

Downstream GWAS workflows use these weights to match GWAS SNPs by chromosome,
position, and allele orientation, then compute PWAS/S-MiXcan association tests.

## Downstream Use

These workflows default to the retained weights file:

```text
Analysis_code/HERMES
Analysis_code/HERMES_HFpEF
Analysis_code/HEARTFAIL
```

The default file used downstream is:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Only change the weights file when intentionally testing a different training
model.

## Useful Test Runs

Run one chromosome only:

```bash
PWAS_CHR_LIST=1 \
PWAS_PRUNE_KB=100 \
PWAS_PRUNE_STEP=1 \
PWAS_PRUNE_R2=0.99 \
bash Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh

PWAS_CHR_LIST=1 \
bash Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh
```

Then combine available chromosomes:

```bash
PWAS_CHR_LIST=1 \
Rscript Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R
```

## Notes

- `1_prune_moderate_heart_protein_dosage.sh` has generic defaults, but the
  retained model uses `100kb`, `step 1`, `r2 = 0.99`.
- The training script does not read GWAS summary statistics.
- GWAS-specific workflows begin after this weight table is available.

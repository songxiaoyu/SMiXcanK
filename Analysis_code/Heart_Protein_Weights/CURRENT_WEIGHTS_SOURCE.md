# Current PWAS Weights Source

The currently retained PWAS training result is:

`Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv`

It was trained by:

1. `Analysis_code/Heart_Protein_Weights/1_prune_moderate_heart_protein_dosage.sh`
2. `Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R`
3. `Analysis_code/Heart_Protein_Weights/run_train_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_by_chr.sh`
4. `Analysis_code/Heart_Protein_Weights/3_combine_heart_protein_moderate_100kb_r2_0.99_alpha0.5_lambdamin_weights.R`

Main training inputs:

- Protein phenotype:
  `Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData`
- Cell fraction estimates:
  `Heart/GTEx_Pi_Estimate/BayesDeBulk_pi.tsv`
- Covariates:
  `New generated files/covariate_EA_with_age.txt`
- EA sample list:
  `Heart/GTEx_Pi_Estimate/GTEx_heart_EA_subject_ids.txt`
- Gene annotation:
  `Data/ensembl38.txt`
- Moderate-pruned GTEx dosage:
  `New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99`

Training settings from the runner:

- Cell components: cardiomyocytes vs other four cell types combined
- Pruning: 100 kb window, `r2 = 0.99`
- MiXcan alpha: `0.5`
- Lambda choice: `lambda.min`
- Output suffix: `_moderate_100kb_r2_0.99_alpha0.5_lambdamin`

HERMES and HERMES_HF prepare scripts now default to this moderate weights file.
Set `PWAS_WEIGHTS_FILE` only when intentionally testing a different weights file.

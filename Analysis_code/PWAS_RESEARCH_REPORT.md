# PWAS / S-MiXcan Research Report

Last updated: 2026-06-29

This document summarizes the analysis history from the original Breast BCAC
two- and three-cell-type S-MiXcan workflows through the newer heart PWAS,
CARDMPRI, HERMES, and diagnostic experiments. It is written as a handoff report
for another agent or analyst to continue debugging why some heart PWAS analyses
produce few or no significant genes.

## 1. Main Question

The central question is:

Why do the heart PWAS analyses produce zero or very few significant genes, while
the older Breast BCAC two- and three-cell-type workflows produce clear signals?

Working hypotheses considered so far:

1. The PWAS training model may have weak predictive power.
2. The PWAS SNP set may be too sparse after LD pruning and overlap filtering.
3. The heart GWAS summary statistics may have weaker signal than BCAC.
4. LD/reference matching or allele matching may reduce usable SNPs.
5. Association-stage regularization may affect stability but is unlikely to
   create true signal if GWAS/model overlap is weak.

## 2. Project Structure

Repository root:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Github
```

Important code folders:

```text
Analysis_code/Breast_BCAC_2CellTypes
Analysis_code/Breast_BCAC_3CellTypes
Analysis_code/Heart_Protein_Weights
Analysis_code/HERMES
Analysis_code/HERMES_HF
```

Important output roots:

```text
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results
/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Figure
```

## 3. Original Breast BCAC Two- and Three-Cell-Type Results

These workflows use breast cell-type training weights and BCAC GWAS summary
statistics. They are the positive-control analyses because they have strong
signals.

### 3.1 Breast BCAC 2-cell-type

Result file:

```text
Results/2pi_workspace/bcac2020_result/bcac2020_result_pi2.csv
```

Summary:

```text
rows tested:        6405
unique genes:       6319
p-value column:     p_join
minimum p:          1.48e-73
Bonferroni < 0.05:  32
FDR < 0.05:         49
FDR < 0.10:         61
median input SNP:   8
mean input SNP:     10.18
genes with 1 SNP:   568
genes with >=5 SNP: 4547
genes with >=10 SNP:2674
```

Top rows:

```text
gene_name          chr    input_snp_num    p_join
ENSG00000182307.12 chr8   4                1.48e-73
ENSG00000179673.4  chr17  2                8.14e-53
ENSG00000072134.15 chr17  4                4.43e-48
ENSG00000139323.13 chr12  11               2.45e-32
```

### 3.2 Breast BCAC 3-cell-type

Result file:

```text
Results/3pi_workspace/bcac2020_result/bcac2020_result_pi3_02.csv
```

Summary:

```text
rows tested:        6403
unique genes:       6320
p-value column:     p_join
minimum p:          1.189e-43
Bonferroni < 0.05:  17
FDR < 0.05:         24
FDR < 0.10:         33
median input SNP:   8
mean input SNP:     10.20
genes with 1 SNP:   669
genes with >=5 SNP: 4481
genes with >=10 SNP:2672
```

Top rows:

```text
gene_name          chr    input_snp_num    p_join
ENSG00000033867.16 chr3   13               1.19e-43
ENSG00000072134.15 chr17  4                1.88e-34
ENSG00000270964.1  chr15  11               3.50e-14
ENSG00000183654.8  chr5   15               3.78e-12
```

Interpretation:

Breast analyses have strong GWAS signal and a moderate SNP count per gene. The
median input SNP count is around 8, and thousands of genes have at least 5 or 10
SNPs.

## 4. PWAS Training Data and Model

The PWAS model uses GTEx heart proteomics with two components:

1. Cardiomyocytes
2. Other, defined as Fibroblasts + Endothelial + Smooth muscle/pericyte + Immune

Main script:

```text
Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R
```

Main protein input:

```text
Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData
```

The protein object uses `gene_id` and `gene_name`; these are Ensembl gene IDs and
gene symbols, not separate protein IDs. Therefore, the PWAS weight table columns
are named `gene_id` and `gene_name`.

Covariate input:

```text
New generated files/covariate_EA_with_age.txt
```

Current training script includes diagnostic output:

```text
cv_mse_cell
cv_mse_tissue
cv_r2_cell
cv_r2_tissue
lambda_cell
lambda_tissue
candidate_snp_num
snp_num_after_maf
selected_snp_num
```

`CV R2` is cross-validation R-squared:

```text
CV R2 = 1 - cross_validated_MSE / variance(centered phenotype)
```

It is used to assess whether genotype weights can predict the protein/cell
component phenotype in held-out folds.

## 5. Original PWAS CARDMPRI Experiment

CARDMPRI GWAS:

```text
Heart/Data/I9_CARDMPRI.gwas.imputed_v3.both_sexes_hg38_rsid.tsv.gz
```

The earlier CARDMPRI setup used strongly pruned genotype dosage:

```text
New generated files/codes/pruned_by_chr
pruning approximately: 500kb, step 1, r2=0.8
```

Main result:

```text
Results/pwas/cardmpri_workspace/cardmpri_result/cardmpri_result_pwas.csv
```

Summary:

```text
rows tested:        5921
unique genes:       5921
p-value column:     p_join
minimum p:          7.023e-04
Bonferroni < 0.05:  0
FDR < 0.05:         0
FDR < 0.10:         0
median input SNP:   1
mean input SNP:     2.72
genes with 1 SNP:   3852
genes with >=5 SNP: 727
genes with >=10 SNP:343
```

Top rows:

```text
gene_name chr   input_snp_num p_join      FDR
STAT1     chr2  1             7.02e-04   0.594
GLS       chr2  2             7.02e-04   0.594
OSGEPL1   chr2  1             7.02e-04   0.594
INPP1     chr2  1             7.02e-04   0.594
```

Interpretation:

This result is highly sparse. The median input SNP per gene is only 1, and most
genes have very few SNPs. This makes the joint S-MiXcan test weak and unstable.

## 6. Archived Original PWAS Experiment

Archived experiment folder:

```text
Results/pwas/experiments/cardmpri_pwas_pruned500kb_r2_0.8_corrected_ncase_20260610
```

This folder contains:

```text
Analysis_code_pwas_snapshot/
figures/
intermediate_outputs/step2_training_model_weights/
intermediate_outputs/step3_cardmpri_input/
intermediate_outputs/step3_4_cardmpri_filtered_id_and_ld/
```

Figures:

```text
figures/CARDMPRI_PWAS_QQ.pdf
figures/CARDMPRI_PWAS_QQ_split_venn.pdf
```

This archive should be treated as the frozen record of the earlier sparse
experiment before changing pruning parameters.

## 7. PWAS Weight Tables

### 7.1 Original sparse weights

File:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other.csv
```

Summary:

```text
rows:               17615
unique genes:       6048
median nonzero SNP: 1
mean nonzero SNP:   2.91
genes with 1 SNP:   3821
genes with >=5 SNP: 783
genes with >=10 SNP:388
```

### 7.2 Moderate pruning weights

To reduce sparsity without using all unpruned SNPs, a moderate pruning strategy
was tested:

```text
100kb window, step 1, r2=0.99
alpha=0.5
lambda choice=lambda.min
```

File:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv
```

Summary:

```text
rows:               154582
unique genes:       6048
median nonzero SNP: 15
mean nonzero SNP:   25.56
genes with 1 SNP:   1573
genes with >=5 SNP: 4049
genes with >=10 SNP:3519
```

Interpretation:

Moderate pruning greatly increases available model SNPs while remaining much
smaller than the fully unpruned genotype matrix.

## 8. Training Diagnostic: chr1 Moderate PWAS

Diagnostic file:

```text
Results/heart_protein_weights/training_model_weights/weights_heart_protein_cardiomyocytes_other_diagnostics_moderate_100kb_r2_0.99_alpha0.5_lambdamin_diagnostic_chr1.csv
```

Run configuration:

```text
PWAS_CHR_FILTER=1
PWAS_GENO_RAW_DIR=New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99
PWAS_MIXCAN_ALPHA=0.5
PWAS_LAMBDA_CHOICE=min
PWAS_WEIGHT_EPS=1e-8
```

Summary:

```text
genes/proteins tested: 632
median CV R2 cell:     0.0917
mean CV R2 cell:       0.1357
median CV R2 tissue:   0.0900
mean CV R2 tissue:     0.1443
cell CV R2 > 0:        490 / 632
cell CV R2 > 0.01:     469 / 632
cell CV R2 > 0.05:     384 / 632
tissue CV R2 > 0:      460 / 632
tissue CV R2 > 0.01:   434 / 632
median selected SNP:   17
median candidate SNP:  2207
median SNP after MAF:  2194.5
selected SNP >=5:      432 / 632
selected SNP >=10:     387 / 632
```

Top chr1 CV R2 examples:

```text
gene_name  selected_snp_num  cv_r2_cell  cv_r2_tissue
CAP1       47                0.737       0.664
GSTM3      29                0.677       0.725
SZRD1      22                0.668       0.676
ACTN2      22                0.604       0.643
HSPB7      33                0.534       0.509
```

Interpretation:

The PWAS training model is not globally failing. On chr1, many proteins have
positive cross-validated R2, and many have nontrivial selected SNP counts. This
supports the idea that the no-signal issue is not simply because the training
model has no predictive ability.

## 9. HERMES DCM PWAS

HERMES input folder:

```text
Heart/HERMES
```

Original HERMES GWAS:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR.tsv.gz
```

Build:

```text
hg19 / GRCh37
```

Case/control numbers:

```text
cases:    14,256
controls: 1,199,156
```

HERMES scripts:

```text
Analysis_code/HERMES/1_liftover_hermes_gwas.py
Analysis_code/HERMES/2_HERMES_prepare_data_pwas.R
Analysis_code/HERMES/3_1000Genome_keep_eur_plink_hermes_pwas.sh
Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
Analysis_code/HERMES/5_run_Primo_hermes_pwas.R
Analysis_code/HERMES/6_plot_hermes_pwas.R
```

The association script includes a tunable inverse regularization parameter:

```text
HERMES_ASSOC_REG_SCALE
default = 0.1
```

This controls regularization during inverse covariance calculation. It may help
numerical stability, but it should not be expected to create biological signal if
SNP overlap or GWAS signal is weak.

### 9.1 HERMES original sparse run

Result file:

```text
Results/hermes_pwas/hermes_workspace/hermes_result/hermes_result_pwas.csv
```

Summary:

```text
rows tested:        5914
unique genes:       5914
p-value column:     p_join
minimum p:          2.272e-05
Bonferroni < 0.05:  0
FDR < 0.05:         0
FDR < 0.10:         0
median input SNP:   1
mean input SNP:     2.50
genes with 1 SNP:   3974
genes with >=5 SNP: 678
genes with >=10 SNP:297
```

Top rows:

```text
gene_name chr   input_snp_num p_join      FDR
SH2D4A    chr8  1             2.27e-05   0.115
PECR      chr2  1             1.12e-04   0.115
IGFBP2    chr2  1             1.12e-04   0.115
IGFBP5    chr2  1             1.12e-04   0.115
```

### 9.2 HERMES moderate pruning run

Workspace:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin
```

Figures:

```text
Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_PWAS_QQ.pdf
Figure/hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_PWAS_QQ_split_venn.pdf
```

Result file:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_result/hermes_result_pwas.csv
```

Summary:

```text
rows tested:        5990
unique genes:       5990
p-value column:     p_join
minimum p:          7.631e-11
Bonferroni < 0.05:  2
FDR < 0.05:         2
FDR < 0.10:         2
median input SNP:   11
mean input SNP:     18.47
genes with 1 SNP:   1578
genes with >=5 SNP: 3902
genes with >=10 SNP:3233
```

Significant genes:

```text
gene_name chr   input_snp_num p_join       FDR
NECAP2    chr1  38            7.63e-11    4.57e-07
HSPB7     chr1  24            2.27e-07    6.79e-04
```

Interpretation:

Using less aggressive pruning increased median input SNP count from 1 to 11 and
produced 2 significant genes in HERMES. This supports the sparsity hypothesis.
However, 2 significant genes is still modest compared with breast BCAC.

## 10. PWAS Experiment 2: PWAS Weights + BCAC GWAS

Purpose:

Use the heart PWAS weights but swap the GWAS to BCAC, to test whether the GWAS
or the PWAS weights are the limiting factor.

Code folder:

```text
Analysis_code/pwas_experiment2
```

Main scripts:

```text
1_BCAC_prepare_data_pwas_weights.R
2_1000Genome_keep_eur_plink_bcac_pwas_weights.sh
3_BCAC_run_analysis_pwas_weights.R
4_run_Primo_bcac_pwas_weights.R
5_plot_bcac_pwas_weights.R
6_diagnostic_compare_gwas_pwas.R
run_step1_to_step5.sh
```

Workspace:

```text
Results/pwas_experiment2/bcac_pwas_weights_workspace
```

Figures:

```text
Figure/pwas_experiment2
```

Result file:

```text
Results/pwas_experiment2/bcac_pwas_weights_workspace/bcac_result/bcac_result_pwas_weights.csv
```

Summary:

```text
rows tested:        5940
unique genes:       5940
p-value column:     p_join
minimum p:          9.645e-07
Bonferroni < 0.05:  1
FDR < 0.05:         1
FDR < 0.10:         1
median input SNP:   1
mean input SNP:     2.55
genes with 1 SNP:   3968
genes with >=5 SNP: 690
genes with >=10 SNP:307
```

Top row:

```text
gene_name chr    input_snp_num p_join      FDR
HTRA1     chr10  18            9.64e-07   0.00573
```

Interpretation:

With the original sparse PWAS weights, BCAC still gives one significant gene, but
far fewer signals than the original Breast BCAC two- and three-cell-type models. This
suggests both factors matter:

1. BCAC is a stronger GWAS than heart GWAS.
2. The heart PWAS weights/overlap are sparse under the original pruning setup.
3. Breast-trained models are better matched to BCAC breast biology than heart
   protein weights are.

## 11. Diagnostic Comparison: Why Breast Has More Signal

Key comparison:

```text
Analysis                  median SNP   FDR<0.05   Bonf<0.05
Breast BCAC 2-cell-type           8            49         32
Breast BCAC 3-cell-type           8            24         17
PWAS CARDMPRI sparse      1            0          0
HERMES sparse             1            0          0
PWAS weights + BCAC       1            1          1
HERMES moderate pruning   11           2          2
```

Important point:

The breast models have a median input SNP count around 8 and strong disease
matched BCAC signal. The original PWAS heart analyses have median input SNP count
1, which is much more sparse. When HERMES is rerun using moderate pruning, median
input SNP increases to 11 and significant genes appear.

Therefore, the strongest current explanation is:

The original PWAS heart analysis was too sparse after SNP pruning and GWAS/LD
overlap filtering. The training model itself has some predictive power, but the
association-stage test loses power when most genes have only 1 or very few SNPs.

## 12. Current Recommendation

Recommended next analysis direction:

1. Continue with moderate pruning rather than the original 500kb r2=0.8 pruning.
   The current reasonable setting is:

```text
100kb window, step 1, r2=0.99
alpha=0.5
lambda=lambda.min
```

2. Re-run CARDMPRI using the moderate PWAS weights and moderate LD/reference
   workspace, then compare against HERMES moderate.

3. Generate a full diagnostic table for CARDMPRI, HERMES, BCAC-control, and
   original BCAC:

```text
analysis
gene/protein count
input_snp_num median
input_snp_num mean
genes with 1 SNP
genes with >=5 SNP
genes with >=10 SNP
min p_join
Bonferroni significant count
FDR < 0.05 count
FDR < 0.10 count
top genes
```

4. For the training model, run CV R2 diagnostics across chr1-22, not only chr1.
   Summarize:

```text
median cv_r2_cell
median cv_r2_tissue
number of proteins with cv_r2_cell > 0
number with cv_r2_cell > 0.01
number with cv_r2_cell > 0.05
relationship between CV R2 and association p-value
```

5. Use association regularization as a stability sensitivity analysis, not as the
   main fix. Suggested values:

```text
HERMES_ASSOC_REG_SCALE=0.01
HERMES_ASSOC_REG_SCALE=0.05
HERMES_ASSOC_REG_SCALE=0.1
HERMES_ASSOC_REG_SCALE=0.2
```

Compare QQ inflation, top genes, and FDR counts across values.

## 13. Open Issues for Next Agent

1. CARDMPRI moderate pruning full run needs to be confirmed or rerun.
2. Full chr1-22 PWAS training diagnostics should be generated.
3. Check whether CARDMPRI case/control numbers are correct. Earlier notes had:

```text
ncase = 360834
ncontrol = 360
```

This looks suspicious because controls are much smaller than cases. It should be
verified against the phenotype definition and summary-stat documentation.

4. Confirm allele orientation and variant ID matching for CARDMPRI and HERMES.
   The scripts attempt forward/reverse allele matching, but this remains an
   important possible failure point.

5. Confirm whether repeated identical minimum p-values in CARDMPRI sparse output
   are expected or indicate an input/association issue.

6. Compare moderate pruning with even less aggressive settings if feasible:

```text
100kb r2=0.995
100kb r2=0.99
50kb r2=0.99
```

The target is to keep the input SNP distribution near or above breast controls
without making matrices too large for the machine.

## 14. Key Takeaway

The evidence so far points away from a complete training-model failure. The
heart PWAS model has measurable CV R2 on chr1, and HERMES begins to show
significant genes once moderate pruning increases SNP coverage. The main current
technical issue is likely sparse SNP overlap after pruning/reference/GWAS
intersection, combined with weaker or phenotype-specific heart GWAS signal.

## 15. Update on 2026-06-29: Package and Workflow Cleanup

The GitHub repository was reorganized so that there is now a single active
top-level R package. The previous nested `SMiXcan/` package copy was moved out
of the active repository and archived externally. The active repository layout
is now:

```text
R/                      package source
man/                    package documentation
data/                   example datasets
Analysis_code/          analysis workflows
  Breast_BCAC_2CellTypes/
  Breast_BCAC_3CellTypes/
  Heart_Protein_Weights/
  HERMES/
  HERMES_HF/
```

The package-level `SMiXcan_assoc_test_K()` was updated. The workflow scripts now
call the package directly with:

```r
library(SMiXcan)
SMiXcan_assoc_test_K(...)
```

They no longer source local copies of `R/S-MiXcan_K.R`.

The active association function supports:

```r
regularization = "fixed"
regularization = "estimate"
```

The old adaptive/tunable association implementation was archived and is not part
of the active workflow. In current documentation, the default package
implementation is treated simply as the association method; no separate "fast"
mode is exposed in the pipeline.

The source package was rebuilt successfully with:

```bash
R CMD build .
```

The final source package excludes `Analysis_code/` through `.Rbuildignore`, so
workflow scripts remain in GitHub but are not included in the R package tarball.

The changes were committed and pushed to GitHub:

```text
commit 1685111
message: Update SMiXcan package and PWAS workflows
remote: git@github.com:songxiaoyu/SMiXcan.git
branch: main
```

This research report was not uploaded because it contains local paths and
internal diagnostic history.

## 16. Final HERMES and HEARTFAIL Results

### 16.1 HERMES DCM PWAS

Input GWAS:

```text
Heart/HERMES/HERMES2_GWAS_DCM_EUR/FORMAT-METAL_Pheno5_EUR.tsv.gz
```

HERMES sample counts:

```text
cases    = 14,256
controls = 1,199,156
```

The case-control effective sample size is:

```text
N_eff = 4 / (1/n_case + 1/n_control) ≈ 56,354
```

Current HERMES result file:

```text
Results/hermes_pwas/hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin_estimate_fast_20260629/hermes_result/hermes_result_pwas.csv
```

Although this workspace name still contains `estimate_fast_20260629`, the active
code now treats the package association implementation as the default method and
does not expose a separate fast mode.

Summary:

```text
rows tested:        5990
usable p-values:    5990
p-value column:     p_join
minimum p:          3.664e-17
median p:           0.477
FDR < 0.10:         5
nominal p < 0.05:   437
median input SNP:   11
mean input SNP:     18.47
genes with 1 SNP:   1578
genes with >=5 SNP: 3902
genes with >=10 SNP:3233
```

Top HERMES rows:

```text
gene_name  chr    input_snp_num  p_join       p_join_pre_reg
NECAP2     chr1   38             3.664e-17    1.336e-06
HSPB7      chr1   24             2.089e-09    3.643e-04
PECR       chr2   15             3.031e-05    4.564e-05
EXOC7      chr17  18             7.190e-05    1.081e-02
PDPK1      chr16  17             7.403e-05    1.598e-04
RPL37A     chr2   1              1.116e-04    1.116e-04
TIMM10     chr11  23             1.552e-04    6.974e-05
RPS27      chr1   22             1.565e-04    6.378e-03
MYBPC3     chr11  44             1.866e-04    3.547e-03
FKBP7      chr2   1              2.420e-04    2.420e-04
```

Interpretation:

HERMES now shows clear evidence of signal. The top results include genes with
cardiac relevance, including `HSPB7` and `MYBPC3`. This supports the view that
the heart protein weights can detect disease-relevant signal when the GWAS has
enough effective sample size and the SNP overlap is not too sparse.

### 16.2 HEARTFAIL PWAS

Input GWAS:

```text
Heart/Data/HEARTFAIL.gwas.imputed_v3.both_sexes.tsv.bgz
```

HEARTFAIL sample counts:

```text
cases    = 1,405
controls = 359,789
```

The case-control effective sample size is:

```text
N_eff = 4 / (1/n_case + 1/n_control) ≈ 5,598
```

Current HEARTFAIL result file:

```text
Results/heartfail_pwas/heartfail_workspace/heartfail_result/heartfail_result_pwas.csv
```

Summary:

```text
rows tested:        5997
usable p-values:    5997
p-value column:     p_join
minimum p:          5.639e-04
median p:           0.504
FDR < 0.10:         0
nominal p < 0.05:   335
median input SNP:   13
mean input SNP:     22.24
genes with 1 SNP:   1554
genes with >=5 SNP: 3987
genes with >=10 SNP:3387
```

Top HEARTFAIL rows:

```text
gene_name    chr    input_snp_num  p_join       p_join_pre_reg
KIF13A       chr6   119            5.639e-04    8.580e-04
MAGOHB       chr12  12             6.307e-04    2.605e-01
VTN          chr17  1              1.035e-03    1.035e-03
SARM1        chr17  1              1.035e-03    1.035e-03
SUPT6H       chr17  1              1.035e-03    1.035e-03
ERAL1        chr17  1              1.035e-03    1.035e-03
TMEM199      chr17  1              1.035e-03    1.035e-03
RP1-66C13.4  chr17  1              1.035e-03    1.035e-03
SDF2         chr17  2              1.116e-03    1.116e-03
ATP6V0A1     chr17  1              1.127e-03    1.127e-03
```

Interpretation:

HEARTFAIL has no genes significant at FDR < 0.10. This does not appear to be
caused by sparse SNP overlap because HEARTFAIL has slightly more input SNPs per
gene than HERMES:

```text
HERMES median input SNP:    11
HEARTFAIL median input SNP: 13
```

The most likely explanation is weaker GWAS power and/or phenotype heterogeneity.
HEARTFAIL has only 1,405 cases and an effective sample size of about 5,598,
roughly one tenth of HERMES. More controls do not compensate for the small number
of cases because case-control GWAS power is limited by the smaller group.

Therefore, HEARTFAIL no-signal is currently interpreted as a GWAS power and
phenotype issue rather than an S-MiXcan implementation failure.

## 17. Effective Sample Size Intercept Simulation

Question tested:

Would scaling the binomial intercept term by case-control effective sample size
improve calibration or recover signal in HERMES?

The relevant intercept term is:

```text
Z0 = log(n_case / n_control) / sqrt(1/n_case + 1/n_control)
```

For HERMES:

```text
Z0 original    ≈ -526.08
N_eff scale    ≈ 0.2155
Z0 effective_n ≈ -113.37
```

A diagnostic script was added locally:

```text
Analysis_code/HERMES/7_simulate_effective_n_intercept.R
```

This script simulates null GWAS Z scores using the real HERMES chr1 weights,
SNP overlap, and LD reference, then compares:

```text
original Z0
effective_n-scaled Z0
```

Simulation run:

```text
chromosome:       chr1
genes tested:     100
simulations:      50
total null tests: 5000
regularization:   estimate
```

Results:

```text
method          n     min p      median p   lambda GC   p<0.05   p<0.01
original        5000  1.779e-05  0.529      0.871       0.0496   0.0106
effective_n     5000  1.779e-05  0.529      0.871       0.0496   0.0106
```

Difference between original and effective_n:

```text
max absolute p difference:    1.25e-12
median absolute p difference: 5.55e-15
max intercept correlation:    6.26e-16
```

Conclusion:

Effective sample size scaling of `Z0` does not materially change the final
SMiXcan p-values. The reason is that the current formula uses:

```r
Y_scaled <- scale(Yhat)
Y <- cbind(1, Y_scaled)
```

After centering `Yhat`, the intercept column is nearly perfectly orthogonal to
the cell-type predicted molecular trait columns. Therefore the large `Z0` from
case-control imbalance does not meaningfully propagate into the cell-type joint
p-values.

This rules out intercept magnitude as the main cause of HEARTFAIL no-signal.

## 18. Current Final Interpretation

The most up-to-date evidence supports the following interpretation:

1. The original heart PWAS analyses were harmed by overly sparse SNP overlap
   after aggressive pruning. Moderate pruning improved HERMES substantially.
2. HERMES has signal after moderate pruning, with 5 genes significant at
   FDR < 0.10.
3. HEARTFAIL still has 0 genes significant at FDR < 0.10 despite slightly better
   SNP coverage than HERMES.
4. The difference between HERMES and HEARTFAIL is most likely driven by GWAS
   power and phenotype specificity:

```text
HERMES effective N:    ~56,354
HEARTFAIL effective N: ~5,598
```

5. Effective sample size scaling of the intercept is not a useful fix because
   it has essentially no effect on the final p-values under the current
   centered design matrix.
6. Further regularization tuning should be treated as sensitivity analysis, not
   as the main explanation for signal/no-signal differences.

Recommended next checks:

```text
1. Compare SNP-level GWAS QQ plots for HERMES and HEARTFAIL.
2. Check whether HEARTFAIL top GWAS SNPs overlap known HF loci.
3. Relate PWAS association p-values to protein model CV R2.
4. Check whether significant HERMES genes concentrate in high-CV-R2 models.
5. If available, rerun HEARTFAIL with a larger or meta-analyzed HF GWAS.
```

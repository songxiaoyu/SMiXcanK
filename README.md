
~~~{r setup, include=FALSE} knitr::opts_chunk\$set( collapse = TRUE,
comment = “\#\>”, message = FALSE, warning = FALSE )~~~

# SMiXcan: Cell-type–aware TWAS for K cell types

**SMiXcan** is a summary statistics-based, cell-type-aware TWAS
framework that infers associations between disease risk and predicted
cell-type-specific expression using only GWAS summary statistics.

This document demonstrates the **basic SMiXcan workflow** in four steps
using **small example data included in the package**.  
The example is intended to illustrate **how to run the functions**, not
to produce biological conclusions.

------------------------------------------------------------------------

## Overview of the workflow

The SMiXcan pipeline consists of four modular steps:

1.  Estimate cell-type fractions from bulk expression  
2.  Train cell-type–specific gene expression prediction models  
3.  Perform SMiXcan association testing  
4.  (Optional) Classify cell-type association patterns using PRIMO

Each step corresponds to one exported function.

------------------------------------------------------------------------

## Step 1: Estimate cell-type fractions

Cell-type fractions represent the proportion of each cell type
contributing to each bulk RNA-seq sample. These fractions can be
estimated from marker genes.

\~\~~{r step1-pi} library(SMiXcanK)

data(exprB_example, markers_example)

res_pi \<- pi_estimation_K( exprB = exprB_example, markers =
markers_example, seed = 1, n.iter = 1000, burn.in = 200 )

head(res_pi\$cell_fraction) \~\~~

### Parameters

- `exprB`  
  Bulk expression matrix (genes × samples).

- `markers`  
  Named list of marker genes defining each cell type. The list length
  determines the number of cell types K.

- `seed`  
  Random seed for reproducibility.

- `n.iter`, `burn.in`  
  MCMC settings passed to the underlying estimator.

### Output

- `cell_fraction`  
  A data frame with one row per sample and one column per cell type.

> SMiXcan does not require this specific estimator.  
> Any N × K matrix of cell-type fractions can be used in later steps.

------------------------------------------------------------------------

## Step 2: Train cell-type–specific MiXcan models

Next, we train genetic prediction models that allow SNP effects on
expression to vary across cell types.

\~\~~{r step2-training} data(x_example, y_example)

pi_k \<- as.matrix(res_pi\$cell_fraction)

fit \<- MiXcan_train_K_symmetric( y = y_example, x = x_example, pi_k =
pi_k, yName = “GENE_A1” )

fit$type fit$beta.SNP.by.cell \~\~~

### Parameters

- `y`  
  Expression vector for a single gene (length = number of samples).

- `x`  
  Genotype matrix for cis-SNPs (samples × SNPs).

- `pi_k`  
  Cell-type fraction matrix (samples × K).

- `yName`  
  Optional gene identifier.

### Output

- `type`  
  Model classification:
  - `"CellTypeSpecific"`
  - `"NonSpecific"`
  - `"NoPredictor"`
- `beta.SNP.by.cell`  
  A list of length K, each element containing SNP weights for one cell
  type.

------------------------------------------------------------------------

## Step 3: SMiXcan association testing

SMiXcan tests gene–trait associations by combining trained SNP weights
with GWAS summary statistics while accounting for LD.

\~\~~{r step3-association} data(gwas_example)

W \<- do.call( cbind,
lapply(fit$beta.SNP.by.cell, function(df) df$weight) )

res_assoc \<- SMiXcan_assoc_test_K( W = W, gwas_results = gwas_example,
x_g = x_example, n0 = 1000, n1 = 1000, family = “binomial” )

res_assoc\$p_join \~\~~

### Parameters

- `W`  
  SNP weight matrix (SNPs × cell types).

- `gwas_results`  
  List with GWAS summary statistics:

  - `Beta`
  - `se_Beta`

- `x_g`  
  Reference genotype matrix used to estimate LD.

- `n0`, `n1`  
  Number of controls and cases.

- `family`  
  Trait type:

  - `"binomial"` (case–control)
  - `"gaussian"` (quantitative)

### Output

- `p_join`  
  Combined association p-value across cell types.

------------------------------------------------------------------------

## Step 4: PRIMO-based pattern classification (optional)

After running SMiXcan genome-wide, results can be post-processed using
PRIMO to classify cell-type association patterns among significant
genes.

\~\~~{r step4-primo} data(merged_example)

res_primo \<- primo_pipeline_wrap( merged = merged_example, pvals_names
= c(“p_1_ct2”, “p_2_ct2”), p_join_name = “p_join_ct2”, type_col =
“type_ct2”, fdr_cutoff = 0.1 )

res_primo$n_shared_total res_primo$tab_specific_patterns \~\~~

### Key outputs

- `out`  
  Input table augmented with:

  - FDR-adjusted p-values
  - PRIMO posterior probabilities
  - MAP pattern labels

- `tab_specific_patterns`  
  Counts of association patterns among significant cell-type-specific
  genes.

- `n_shared_total`  
  Total number of shared genes.

------------------------------------------------------------------------

## Summary

SMiXcan provides a flexible framework to:

- model gene expression prediction across K cell types
- allow shared and cell-type–specific genetic effects
- perform joint TWAS using GWAS summary statistics
- optionally classify cell-type association patterns using PRIMO

Each step is modular and can be replaced or extended as needed.

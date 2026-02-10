
# SMiXcan: Cell-type–aware TWAS for K cell types

**SMiXcan** is a summary statistics-based, cell-type-aware TWAS
framework that infers associations between disease risk and predicted
cell-type-specific expression using only GWAS summary statistics.

This document demonstrates the **basic SMiXcan workflow** in four steps
using **small example data included in the package**.  
The example is intended to illustrate **how to run the functions**, not
to produce biological conclusions.

------------------------------------------------------------------------

``` r
knitr::opts_chunk$set(echo = TRUE)
library(SMiXcanK) 
```

## Overview of the workflow

The SMiXcan pipeline consists of four modular steps:

1.  Estimate cell-type fractions from bulk expression  
2.  Train cell-type–specific gene expression prediction models  
3.  Perform SMiXcan association testing  
4.  (Optional) Classify cell-type association patterns using PRIMO

Each step corresponds to one exported function.

------------------------------------------------------------------------

## Package Installation

With R, users can install the SMiXcan package directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/SMiXcanK")
```

------------------------------------------------------------------------

## Step 1: Estimate cell-type fractions

Cell-type fractions represent the proportion of each cell type
contributing to each bulk RNA-seq sample. These fractions can be
estimated from marker genes.

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

``` r
library(SMiXcanK) 
data(exprB_example, markers_example) 
head(exprB_example)
```

    ##         Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
    ## GENE_A1     8.1     8.0     7.9     8.2     8.1    7.95    8.05    8.15
    ## GENE_A2     7.5     7.6     7.4     7.7     7.6    7.55    7.45    7.65
    ## GENE_A3     6.9     7.0     6.8     7.1     7.0    6.95    6.85    7.05
    ## GENE_B1     5.2     5.1     5.3     5.2     5.1    5.15    5.25    5.05
    ## GENE_B2     4.8     4.7     4.9     4.8     4.7    4.75    4.85    4.65
    ## GENE_B3     4.3     4.4     4.2     4.5     4.4    4.35    4.25    4.45

``` r
head(markers_example)
```

    ## $CellType1
    ## [1] "GENE_A1" "GENE_A2" "GENE_A3"
    ## 
    ## $CellType2
    ## [1] "GENE_B1" "GENE_B2" "GENE_B3"

``` r
res_pi <- pi_estimation_K( exprB = exprB_example, markers = markers_example, seed = 1, n.iter = 1000, burn.in = 200 )
head(res_pi$cell_fraction) 
```

    ##   SampleID  CellType1  CellType2
    ## 1  Sample1 0.28802059 0.71197941
    ## 2  Sample2 0.87463867 0.12536133
    ## 3  Sample3 0.06539837 0.93460163
    ## 4  Sample4 0.89385834 0.10614166
    ## 5  Sample5 0.98156287 0.01843713
    ## 6  Sample6 0.77057284 0.22942716

### Output

- `cell_fraction`  
  A data frame with one row per sample and one column per cell type.

> SMiXcan does not require this specific estimator.  
> Any N × K matrix of cell-type fractions can be used in later steps.

------------------------------------------------------------------------

## Step 2: Train cell-type–specific MiXcan models

Next, we train genetic prediction models that allow SNP effects on
expression to vary across cell types using a **symmetric
parameterization**.

### Parameters

- `y`  
  Expression vector for a single gene (length = number of samples).

- `x`  
  Genotype matrix for cis-SNPs (samples × SNPs).

- `pi_k`  
  Cell-type fraction matrix (samples × K).

- `yName`  
  Optional gene identifier.

``` r
data(x_example, y_example)
head(x_example)
```

    ##         rs1 rs2 rs3 rs4 rs5
    ## Sample1   0   1   2   0   1
    ## Sample2   1   0   1   2   0
    ## Sample3   2   1   0   1   2
    ## Sample4   0   2   1   0   1
    ## Sample5   1   1   2   1   0
    ## Sample6   2   0   1   2   1

``` r
head(y_example)
```

    ## Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 
    ##    2.10    2.30    1.90    2.20    2.00    2.15

``` r
pi_k <- as.matrix(res_pi$cell_fraction)
rownames(pi_k) <- pi_k[, 1]          
pi_k <- pi_k[, -1, drop = FALSE]    
pi_k <- apply(pi_k, 2, as.numeric)   
pi_k <- as.matrix(pi_k)  

foldid <- sample(rep(1:4, length.out = nrow(x_example)))
fit <- MiXcan_train_K_symmetric(
  y     = y_example,
  x     = x_example,
  pi_k  = pi_k,
  foldid = foldid,
  yName = "GENE_A1"
)
```

    ## Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per
    ## fold
    ## Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per
    ## fold

``` r
fit$type 
```

    ## [1] "NonSpecific"

``` r
fit$beta.SNP.by.cell
```

    ## $Cell1
    ##     SNP     weight
    ## rs1 rs1 -0.1106899
    ## rs2 rs2  0.0000000
    ## rs3 rs3  0.0000000
    ## rs4 rs4  0.0000000
    ## rs5 rs5  0.0000000
    ## 
    ## $Cell2
    ##     SNP     weight
    ## rs1 rs1 -0.1106899
    ## rs2 rs2  0.0000000
    ## rs3 rs3  0.0000000
    ## rs4 rs4  0.0000000
    ## rs5 rs5  0.0000000

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

``` r
data(gwas_example)

W <- do.call( cbind, lapply(fit$beta.SNP.by.cell, function(df) df$weight) )

res_assoc <- SMiXcan_assoc_test_K( W = W, gwas_results = gwas_example, x_g = x_example, n0 = 1000, n1 = 1000, family = "binomial" )

res_assoc$p_join 
```

    ## [1] 5.733031e-07

### Output

- `p_join`  
  Combined association p-value across cell types.

## Summary

SMiXcan provides a flexible framework to:

- model gene expression prediction across K cell types
- perform cell-type-aware TWAS using GWAS summary statistics
- optionally classify cell-type association patterns using PRIMO

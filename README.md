
# SMiXcan: Cell-type–aware TWAS from Bulk Transcriptomics Using GWAS Summary Statistics for K Cell Types

**SMiXcan** is a summary statistics-based, cell-type-aware TWAS
framework that infers associations between disease risk and predicted
cell-type-specific expression using only GWAS summary statistics.

This document demonstrates the **basic SMiXcan workflow** in four steps
using **small example data included in the package**. The example is
intended to illustrate **how to run the functions**, not to produce
biological conclusions.

------------------------------------------------------------------------

## Overview of the workflow

The SMiXcan pipeline consists of four modular steps:

1.  Estimate cell-type fractions from bulk expression
2.  Train cell-type–specific gene expression prediction models
3.  Perform SMiXcan association testing
4.  (Optional) Classify cell-type association patterns

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

## Step 1: Estimate cell-type fractions (Optional)

Cell-type fractions represent the proportion of each cell type
contributing to each bulk RNA-seq sample. These fractions can be
estimated from marker genes. This step is optional. Users can use their
own cell-type proportion estimates instead.

### Parameters

- `exprB` Bulk expression matrix (genes × samples).

- `markers` Named list of marker genes defining each cell type. The list
  length determines the number of cell types K.

- `seed` Random seed for reproducibility.

- `n.iter`, `burn.in` MCMC settings passed to the underlying estimator.

``` r
library(SMiXcanK)
data(exprB_example, markers_example)
head(exprB_example)
```

    ##          Sample1    Sample2    Sample3    Sample4     Sample5    Sample6
    ## Gene1 -0.6357365  0.4251004 -0.1643758  0.9510128  1.04361246 -2.5923277
    ## Gene2 -0.4616447 -0.2386471  0.4206946 -0.3892372  0.09907849  1.3140022
    ## Gene3  1.4322822  1.0584830 -0.4002467 -0.2843307 -0.45413691 -0.6355430
    ## Gene4 -0.6506964  0.8864227 -1.3702079  0.8574098 -0.65578185 -0.4299788
    ## Gene5 -0.2073807 -0.6192430  0.9878383  1.7196273 -0.03592242 -0.1693183
    ## Gene6 -0.3928079  2.2061025  1.5197450  0.2700549  1.06916146  0.6122182
    ##          Sample7     Sample8     Sample9    Sample10   Sample11   Sample12
    ## Gene1 -0.3881675 -0.45303708 -2.52850069 -1.21536404  0.1643729  0.3309763
    ## Gene2  0.6525365  2.16536850 -0.93590256 -0.02255863  0.3326236  0.9763275
    ## Gene3  1.1247724  1.24574667 -0.96723946  0.70123930 -0.3852080 -0.8433399
    ## Gene4 -0.7721108  0.59549803  0.04748859 -0.58748203 -1.3987540 -0.9705799
    ## Gene5 -0.5080862  0.00488445 -0.40373679 -0.60672794  2.6757408 -1.7715313
    ## Gene6  0.5236206  0.27936078  0.23149613  1.09664022 -0.4236861 -0.3224703
    ##         Sample13    Sample14   Sample15   Sample16    Sample17    Sample18
    ## Gene1  1.5683647 -0.01128123  0.3533985  0.7960927  0.09535401 -0.79506319
    ## Gene2  1.2967556  0.61967726  0.7268137  0.9864283 -0.46281942 -0.01995512
    ## Gene3 -0.2375963 -1.28123874  0.6682610 -0.7945317 -1.46888216 -2.51442512
    ## Gene4 -1.2241501 -0.12426133 -2.4243173 -0.3088180  0.15268651  2.21095203
    ## Gene5 -0.3278127  0.17574165 -0.2353574  0.3614448  1.77376261 -1.48876223
    ## Gene6 -2.4124503  1.69277379  1.9796333  1.3987911 -0.64807093 -1.16075188
    ##         Sample19   Sample20
    ## Gene1 -0.5408727 -2.3214909
    ## Gene2 -0.2163758  1.3641192
    ## Gene3 -1.6219373  1.1322291
    ## Gene4 -1.4509640 -0.7743163
    ## Gene5  0.3509097 -1.4103750
    ## Gene6 -0.1745469 -1.8345276

``` r
head(markers_example)
```

    ## $CellType1
    ## [1] "Gene1" "Gene2" "Gene3"
    ## 
    ## $CellType2
    ## [1] "Gene4" "Gene5" "Gene6"

``` r
res_pi <- pi_estimation_K(exprB = exprB_example, markers = markers_example, seed = 1, n.iter = 1000, burn.in = 200)
head(res_pi$cell_fraction)
```

    ##         CellType1 CellType2
    ## Sample1 0.8539589 0.1460411
    ## Sample2 0.1989060 0.8010940
    ## Sample3 0.7457105 0.2542895
    ## Sample4 0.1908111 0.8091889
    ## Sample5 0.5677247 0.4322753
    ## Sample6 0.3247569 0.6752431

### Output

- `cell_fraction` A data frame with one row per sample and one column
  per cell type.

------------------------------------------------------------------------

## Step 2: Train cell-type–specific MiXcan models

Next, we train genetic prediction models that allow SNP effects on
expression to vary across cell types using a symmetric parameterization.

### Parameters

- `y` Expression vector for a single gene (length = number of samples).

- `x` Genotype matrix for cis-SNPs (samples × SNPs).

- `pi_k` Cell-type fraction matrix (samples × K).

- `yName` Optional gene identifier.

``` r
data(x_example)
data(y_example)
data(pi_k)

set.seed(1)
foldid <- sample(rep(1:4, length.out = nrow(x_example)))

fit <- MiXcan_train_K(
  y      = y_example,
  x      = x_example,
  pi_k   = pi_k,
  foldid = foldid,
  yName  = "GENE_A1"
)

fit$type
```

    ## [1] "CellTypeSpecific"

``` r
fit$W
```

    ##               Cell1         Cell2
    ## SNP1  -4.192148e-01 -4.192148e-01
    ## SNP2   7.557131e-17 -7.557131e-17
    ## SNP3   0.000000e+00  0.000000e+00
    ## SNP4   0.000000e+00  0.000000e+00
    ## SNP5   0.000000e+00  0.000000e+00
    ## SNP6   0.000000e+00  0.000000e+00
    ## SNP7   0.000000e+00  0.000000e+00
    ## SNP8   0.000000e+00  0.000000e+00
    ## SNP9   0.000000e+00  0.000000e+00
    ## SNP10  0.000000e+00  0.000000e+00

### Output of `MiXcan_train_K()`

The object `fit` returned by `MiXcan_train_K()` is a list containing:

- `type`

  - “CellTypeSpecific”
  - “NonSpecific”
  - “NoPredictor”

- `beta.SNP.by.cell`

- `beta.all.models`

- `W`

- `glmnet.cell`

- `glmnet.tissue`

- `yName`

- `xNameMatrix`

------------------------------------------------------------------------

## Step 3: SMiXcan association testing

SMiXcan tests gene–trait associations by combining trained SNP weights
with GWAS summary statistics while accounting for LD.

### Parameters

- `W`
- `gwas_results`
- `x_g`
- `n0`, `n1`
- `family`

``` r
W <- fit$W

res_assoc <- SMiXcan_assoc_test_K(
  W            = W,
  gwas_results = gwas_example,
  x_g          = x_example,
  n0           = 1000,
  n1           = 1000,
  family       = "binomial"
)

res_assoc
```

    ## $Z_join
    ## [1] -3.084652 -3.084652
    ## 
    ## $p_join_vec
    ## [1] 0.002037907 0.002037907
    ## 
    ## $p_join
    ## [1] 0.002037907

### Output of `SMiXcan_assoc_test_K()`

- `Z_join`
- `p_join_vec`
- `p_join`

### Interpretation

- If `p_join` is significant, the gene shows evidence of association.
- `p_join_vec` indicates which cell types contribute most strongly.
- The sign of `Z_join` indicates direction of effect.

------------------------------------------------------------------------

## Step 4: PRIMO-based cell-type association pattern classification (optional)

After running SMiXcan genome-wide, we may wish to classify significant
genes into cell-type association patterns using a PRIMO-based framework.

This step:

1.  Adjusts marginal and joint p-values using BH-FDR  
2.  Identifies significant genes  
3.  Estimates PRIMO posterior probabilities  
4.  Assigns MAP patterns among non-null configurations

``` r
library(SMiXcanK)

# Load example merged results
data(merged_example)

head(merged_example)
```

    ##    gene         type_ct2      p_1_ct2 p_2_ct2 p_join_ct2
    ## 1 Gene1 CellTypeSpecific 1.000000e-08     0.5        0.5
    ## 2 Gene2 CellTypeSpecific 1.623777e-08     0.5        0.5
    ## 3 Gene3 CellTypeSpecific 2.636651e-08     0.5        0.5
    ## 4 Gene4 CellTypeSpecific 4.281332e-08     0.5        0.5
    ## 5 Gene5 CellTypeSpecific 6.951928e-08     0.5        0.5
    ## 6 Gene6 CellTypeSpecific 1.128838e-07     0.5        0.5

``` r
res_primo <- primo_pipeline_wrap(
  merged          = merged_example,
  pvals_names     = c("p_1_ct2", "p_2_ct2"),
  p_join_name     = "p_join_ct2",
  type_col        = "type_ct2",
  fdr_cutoff      = 0.1
)
```

    ## 
    ## Iteration: 1
    ## Iteration: 2
    ## Iteration: 3

``` r
res_primo$n_trait
```

    ## cell1 cell2 
    ##    22    20

``` r
res_primo$n_shared_specific
```

    ## [1] 18

``` r
res_primo$n_shared_nonspecific
```

    ## [1] 15

``` r
res_primo$n_shared_total
```

    ## [1] 33

``` r
res_primo$tab_specific_patterns
```

    ## 
    ## 01 10 11 
    ## 20 22 18

### Parameters

- `merged`  
  Data frame containing:

  - marginal p-values for each cell type  
  - joint p-value  
  - cell-type specificity label

- `pvals_names`  
  Character vector of marginal p-value column names (length = K).

- `p_join_name`  
  Column name of joint (non-specific) p-value.

- `type_col`  
  Column indicating `"CellTypeSpecific"` vs `"NonSpecific"` genes.

- `fdr_cutoff`  
  FDR threshold used to define significance (default 0.1).

------------------------------------------------------------------------

### Key Outputs

- `out`  
  Input table augmented with:

  - BH-adjusted p-values  
  - PRIMO posterior probabilities (`post_*`)  
  - MAP pattern labels (`MAP_pattern_nonnull`)

- `n_trait`  
  Number of genes uniquely associated with each cell type.

- `n_shared_specific`  
  Number of significant cell-type-specific genes associated with
  multiple cell types.

- `n_shared_nonspecific`  
  Number of significant non-specific genes (based on joint test).

- `n_shared_total`  
  Total number of shared genes.

- `tab_specific_patterns`  
  Table summarizing MAP pattern counts among significant
  cell-type-specific genes.

### Interpretation

For K = 2 cell types, PRIMO considers the following patterns:

- `"10"` : associated only with cell type 1  
- `"01"` : associated only with cell type 2  
- `"11"` : associated with both cell types

MAP classification is performed **only among significant
cell-type-specific genes** and excludes the null pattern (`"00"`).

------------------------------------------------------------------------

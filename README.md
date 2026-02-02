
# S-MiXcan

**S-MiXcan** is an R package for performing cell-type-aware
**summary-based transcriptome-wide association studies (TWAS)**. It
extends the MiXcan framework to analyze associations between genetically
regulated gene expression (GReX) and traits and their functional cell
types using GWAS summary statistics.

------------------------------------------------------------------------

## 🔬 Overview

Traditional TWAS approaches predict gene expression at the *bulk tissue*
level and test its association with disease traits, ignoring
heterogeneity across cell types. **MiXcan** improves on this by enabling
cell-type-aware TWAS, but it requires *individual-level genotype data*.

**S-MiXcan** addresses this limitation.  
It provides a **summary-statistics-based** framework that:

- Infers cell-type-specific GReX–trait associations from GWAS summary
  statistics.
- Leverages MiXcan-trained cell-type models.
- Evaluates functional cell-types.
- Requires **no individual-level genotype data**.

------------------------------------------------------------------------

## 📦 Installation

You can install the development version of `S-MiXcan` from GitHub:

``` r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install S-MiXcan from GitHub
devtools::install_github("songxiaoyu/SMiXcanK")

library(SMiXcanK)
```

------------------------------------------------------------------------

## 🧪 Example Usage

``` r
# Define input 
## Example: regularized_inverse_cov()

set.seed(1)

# Create a well-behaved covariance matrix (p x p)
p <- 5
A <- matrix(rnorm(p * p), p, p)
X <- crossprod(A)  # symmetric positive-definite covariance-like matrix

# Regularized inverse
ri <- regularized_inverse_cov(X)

# Outputs
ri$lambda        # regularization strength used (scalar here)
```

    ##            [,1]       [,2]       [,3]       [,4]       [,5]
    ## [1,] 0.10000000 0.06833590 0.08766259 0.06451824 0.09455443
    ## [2,] 0.06833590 0.10000000 0.09423659 0.08650160 0.08679522
    ## [3,] 0.08766259 0.09423659 0.10000000 0.08587295 0.09832630
    ## [4,] 0.06451824 0.08650160 0.08587295 0.10000000 0.08199720
    ## [5,] 0.09455443 0.08679522 0.09832630 0.08199720 0.10000000

``` r
ri$inv[1:3, 1:3] # top-left corner of the inverse matrix
```

    ##             [,1]       [,2]        [,3]
    ## [1,]  0.44205757 0.02250073 -0.01685235
    ## [2,]  0.02250073 1.19681386  0.54483091
    ## [3,] -0.01685235 0.54483091  0.63059602

``` r
## Example: SMiXcan_assoc_test_K() with K = 3 cell types

set.seed(123)

# Dimensions
n  <- 400   # reference panel individuals
p  <- 60    # SNPs in gene region
K  <- 3     # cell types

# Reference genotype matrix (n x p); here simulated as standardized normals
x_g <- scale(matrix(rnorm(n * p), n, p))

# Cell-type weights matrix (p x K)
W <- matrix(rnorm(p * K), p, K)

# GWAS summary stats for the same p SNPs
gwas_results <- list(
  Beta    = rnorm(p, mean = 0, sd = 0.05),
  se_Beta = runif(p, min = 0.02, max = 0.08)
)

# Case/control sample sizes (for binomial null model)
n1 <- 5000
n0 <- 5000

# Run S-MiXcan association test (K = 3)
res <- SMiXcan_assoc_test_K(
  W = W,
  gwas_results = gwas_results,
  x_g = x_g,
  n0 = n0,
  n1 = n1,
  family = "binomial"
)

# Results:
# - Z_join: joint-model Z for each cell type (length K)
# - p_join_vec: per-cell-type joint p-values (length K)
# - p_join: ACAT-combined p across K cell types
res$Z_join
```

    ## [1] 0.06863635 1.04781404 0.69442659

``` r
res$p_join_vec
```

    ## [1] 0.9452791 0.2947243 0.4874147

``` r
res$p_join
```

    ## [1] 0.827071

------------------------------------------------------------------------

## 📄 Citation

If you use **S-MiXcan** in your research, please cite this page

------------------------------------------------------------------------

## 📫 Contact

For questions, please contact:

**Sinan Zhu**  
PhD Candidate, Duke-NUS Medical School  
Email: <sinan.zhu@u.duke.nus.edu>

------------------------------------------------------------------------

## 🔒 License

------------------------------------------------------------------------

\`\`\`r rmarkdown::render(“README.Rmd”)

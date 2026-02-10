#' Example bulk expression matrix
#'
#' A small toy bulk expression matrix (genes x samples) used in the README.
#'
#' @format A numeric matrix with genes in rows and samples in columns.
#' @source Included with the SMiXcanK package.
"exprB_example"

#' Example marker gene list
#'
#' Marker genes for each cell type used by `pi_estimation_K()`.
#'
#' @format A named list of character vectors.
#' @source Included with the SMiXcanK package.
"markers_example"

#' Example genotype matrix
#'
#' A small toy genotype matrix (samples x SNPs) used in the README.
#'
#' @format A numeric matrix.
#' @source Included with the SMiXcanK package.
"x_example"

#' Example gene expression vector
#'
#' A small toy expression vector for one gene (length = samples).
#'
#' @format A numeric vector.
#' @source Included with the SMiXcanK package.
"y_example"

#' Example GWAS summary statistics
#'
#' A minimal list containing GWAS effect sizes and standard errors for the SNPs
#' in `x_example`.
#'
#' @format A list with elements `Beta` and `se_Beta`.
#' @source Included with the SMiXcanK package.
"gwas_example"

#' Example merged table for PRIMO post-processing
#'
#' A toy results table with marginal p-values per cell type, a joint p-value, and
#' a specificity label, used to demonstrate `primo_pipeline_wrap()`.
#'
#' @format A data.frame.
#' @source Included with the SMiXcanK package.
"merged_example"

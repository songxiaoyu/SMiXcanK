#' Example data for S-MiXcan
#'
#' These datasets are provided as examples for demonstrating the usage of
#' the S-MiXcan package. They include example genotype reference data,
#' GWAS summary statistics, and phenotype distributions.
#'
#' @format Depends on the dataset:
#'
#' \describe{
#'   \item{W}{A matrix of cell-type-specific weights (example).}
#'   \item{X_ref_filtered}{Reference genotype matrix after LD filtering.}
#'   \item{gwas_results}{Example GWAS summary statistics.}
#'   \item{n0}{Number of controls for GWAS sample size (example).}
#'   \item{n1}{Number of cases for GWAS sample size (example).}
#' }
#'
#' @source Generated for demonstration purposes.
#'
"W"

#' Example reference genotype matrix
#'
#' An example reference genotype matrix after quality control and
#' LD pruning, used for demonstrating S-MiXcan.
#'
#' @format A numeric matrix with SNPs in columns and individuals in rows.
#' @docType data
#' @usage data(X_ref_filtered)
#' @keywords datasets
"X_ref_filtered"


#' Example GWAS summary statistics
#'
#' Example GWAS summary statistics used to demonstrate S-MiXcan
#' association testing.
#'
#' @format A data frame with SNP-level summary statistics, including
#' columns such as SNP ID, effect allele, non-effect allele, effect
#' size estimate, standard error, and p-value.
#' @docType data
#' @usage data(gwas_results)
#' @keywords datasets
"gwas_results"


#' Example number of controls
#'
#' Example value for the number of controls in a case-control GWAS.
#'
#' @format A single numeric value.
#' @docType data
#' @usage data(n0)
#' @keywords datasets
"n0"


#' Example number of cases
#'
#' Example value for the number of cases in a case-control GWAS.
#'
#' @format A single numeric value.
#' @docType data
#' @usage data(n1)
#' @keywords datasets
"n1"


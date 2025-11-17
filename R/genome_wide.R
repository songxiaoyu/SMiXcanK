#' Run GWAS for Each SNP
#'
#' Performs univariate GWAS by regressing the phenotype on each SNP separately.
#'
#' @name run_gwas
#' @title Run Genome-Wide Association Study
#'
#' @param X A genotype matrix (individuals × SNPs).
#' @param D A phenotype vector (numeric or binary).
#' @param family0 A GLM family, either "binomial" or "gaussian".
#'
#' @return A list with:
#' \describe{
#'   \item{Beta}{Vector of effect sizes for each SNP}
#'   \item{se_Beta}{Vector of standard errors for each SNP}
#' }
#'
#' @importFrom tibble lst
#' @export

run_gwas <- function(X, D, family0) {
  p <- ncol(X)
  Beta <- numeric(p)
  se_Beta <- numeric(p)

  for (j in seq_len(p)) {
    fit <- stats::glm(D ~ X[, j], family = family0)
    cs <- summary(fit)$coefficients
    Beta[j]   <- cs[2, "Estimate"]
    se_Beta[j] <- cs[2, "Std. Error"]
  }

  results <- list(Beta = Beta, se_Beta = se_Beta)
  return(results)
}

#' Run GWAS in parallel across SNPs or predictors
#'
#' @description
#' Fits univariate GLMs for each column of X in parallel using foreach/doParallel,
#' returning the estimated effect size (beta) and standard error for each SNP or predictor.
#'
#' @param X A numeric matrix or data frame of predictors (n × p).
#' @param D A numeric vector or factor of outcomes (length n).
#' @param family0 A valid family object (e.g. \code{gaussian()}, \code{binomial()}).
#'
#' @return A list with:
#' \item{Beta}{Estimated regression coefficients for each column of X.}
#' \item{se_Beta}{Standard errors for each coefficient.}
#'
#' @importFrom stats glm
#' @export
run_gwas_parallel <- function(X, D, family0) {
  r1 <- foreach::foreach(j = 1:ncol(X), .combine = "rbind") %dopar% {
    gwas_result <- stats::glm(D ~ X[, j], family = family0)
    coef_summary <- summary(gwas_result)$coefficients
    coef_summary[2, 1:2]  # beta and SE
  }

  results <- list(
    Beta = r1[, 1],
    se_Beta = r1[, 2]
  )
  return(results)
}


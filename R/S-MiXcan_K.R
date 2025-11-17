#' Regularized Inverse of a Covariance Matrix
#'
#' Computes a regularized inverse of a covariance matrix using correlation shrinkage.
#'
#' @name regularized_inverse_cov
#' @title Regularized Covariance Matrix Inversion
#'
#' @param X A numeric covariance matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{inv}{The regularized inverse matrix.}
#'   \item{lambda}{The regularization parameter used.}
#' }
#'
#' @export

regularized_inverse_cov <- function(X) {

  r <- abs(cor(X))
  lambda = 0.001 * abs(r)^6
  # Apply regularization if lambda > 0
  X_cor <- cov2cor(X)
  X_reg <- X_cor + lambda * diag(nrow(X))
  X_cor_inv <- solve(X_reg)
  D <- diag(X)
  X_inv <- diag(1/sqrt(D)) %*% X_cor_inv %*% diag(1/sqrt(D))

  return(list(
    inv = X_inv,
    lambda = lambda
  ))
}

#' S-MiXcan Association Test with Shrinkage
#'
#' Performs the S-MiXcan association test using GWAS summary statistics and reference genotype data,
#' applying shrinkage-based regularization to stabilize inverse covariance estimation.
#'
#' @param W A p-by-k matrix of cell-type level weights  (where p is the number of SNPs in the gene region, k is the number of cell types).
#' @param gwas_results A list containing GWAS summary statistics, with components \code{Beta} and \code{se_Beta}.
#' @param x_g A genotype matrix for the reference panel (individuals × SNPs).
#' @param n0 Number of controls.
#' @param n1 Number of cases.
#' @param family Either \code{"binomial"} or \code{"gaussian"} (used for fitting the null model).
#' @return A list containing:
#' \describe{
#'   \item{Z_join}{Z-score for cell-type 1 to K (joint model)}
#'   \item{p_join_vec}{P-value for cell-type 1 to K (joint model)}
#'   \item{p_join}{Combined p-value from joint model (ACAT)}
#' }
#'
#' @importFrom tibble lst
#' @export
SMiXcan_assoc_test_K <- function(W,
                                 gwas_results,
                                 x_g,
                                 n0,
                                 n1,
                                 family = c("binomial", "gaussian")) {
  family <- match.arg(family)

  # ---- Inputs and basic checks ----
  Beta    <- as.numeric(gwas_results$Beta)
  se_Beta <- as.numeric(gwas_results$se_Beta)

  if (!is.matrix(W))
    stop("W must be a numeric matrix of dimension p × K.")

  p <- nrow(W)
  K <- ncol(W)

  if (ncol(x_g) != p)
    stop("Number of SNPs (columns) in x_g must match nrow(W).")

  if (length(Beta) != p || length(se_Beta) != p)
    stop("Length of Beta and se_Beta must match nrow(W).")

  # ---- LD and SNP variances ----
  cov_x <- stats::var(x_g)           # p × p LD (covariance) matrix
  sig_l <- sqrt(diag(cov_x))         # per-SNP SD

  # ---- Separate-component Z for each cell type ----
  Z_sep <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    wk <- W[, k]
    sig2_gk <- drop(t(wk) %*% cov_x %*% wk)  # Var(predicted expression)
    if (sig2_gk > 0) {
      num <- sum((wk * Beta) * sig_l / se_Beta)
      Z_sep[k] <- num / sqrt(sig2_gk)
    } else {
      Z_sep[k] <- NA_real_
    }
  }

  p_sep <- ifelse(
    is.na(Z_sep),
    NA_real_,
    2 * stats::pnorm(abs(Z_sep), lower.tail = FALSE)
  )

  # Default joint outputs = separate
  Z_join     <- Z_sep
  p_join_vec <- p_sep
  mode       <- "separate"

  # ---- Joint test(s) ----
  if (family == "binomial") {
    # 1. Null intercept-only logistic => Z0
    D   <- c(rep(1L, n1), rep(0L, n0))
    fit0 <- stats::glm(D ~ 1, family = stats::binomial())
    coef0 <- summary(fit0)$coefficients
    Z0    <- coef0["(Intercept)", "z value"]

    # 2. Predicted expression for each cell type
    Yhat <- x_g %*% W             # n × K
    Y_scaled <- scale(Yhat)       # standardize columns
    Y <- cbind(1, Y_scaled)       # [Intercept, cell-type 1..K]
    colnames(Y) <- c("Y0", paste0("Y", seq_len(K)))

    YtY   <- crossprod(Y)         # (K+1) × (K+1)
    Omega <- diag(YtY)

    # Correlation matrix on the Y scale
    corY <- stats::cov2cor(YtY)

    # Near-singularity check: any |cor| > threshold?
    if (!any(is.na(corY)) &&
        !any(abs(corY[upper.tri(corY)]) > 0.999999)) {

      S <- regularized_inverse_cov(YtY)$inv
      v <- c(
        sqrt(Omega[1]) * Z0,
        sqrt(Omega[-1]) * Z_sep
      )

      Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
      Z_join <- as.numeric(Z_full[-1])  # drop intercept
      p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
      mode <- "joint"
    }

  } else if (family == "gaussian") {
    # Gaussian joint step uses only cell-type components (no intercept)
    Yhat <- x_g %*% W
    Y_scaled <- scale(Yhat)             # n × K
    colnames(Y_scaled) <- paste0("Y", seq_len(K))

    YtY   <- crossprod(Y_scaled)        # K × K
    Omega <- diag(YtY)
    corY  <- stats::cov2cor(YtY)

    if (!any(is.na(corY)) &&
        !any(abs(corY[upper.tri(corY)]) > 0.999999)) {

      S <- regularized_inverse_cov(YtY)$inv
      v <- sqrt(Omega) * Z_sep          # length K

      Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
      Z_join <- as.numeric(Z_full)
      p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
      mode <- "joint"
    }

  } else {
    stop("family must be 'binomial' or 'gaussian'")
  }

  # ---- ACAT-combined p-value across K cell types ----
  p_join <- safe_ACAT(p_join_vec)

  list(
    Z_join    = Z_join,
    p_join_vec = p_join_vec,
    p_join    = p_join,
  )
}



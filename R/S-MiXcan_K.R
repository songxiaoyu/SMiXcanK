#' Compute a regularized inverse covariance matrix
#'
#' This helper stabilizes covariance-matrix inversion by applying diagonal
#' regularization to the correlation matrix before inversion.
#'
#' @name regularized_inverse_cov
#' @title Regularized Covariance Matrix Inversion
#'
#' @param X Numeric covariance matrix.
#' @param reg_scale Non-negative diagonal regularization scale.
#'
#' @return A list containing:
#' \describe{
#'   \item{inv}{The regularized inverse matrix.}
#'   \item{lambda}{The regularization parameter used.}
#'   \item{condition}{Condition number after regularization.}
#' }
#'
#' @export
regularized_inverse_cov <- function(X, reg_scale = 0.1) {
  X_cor <- stats::cov2cor(X)
  X_reg <- X_cor + reg_scale * diag(nrow(X_cor))
  X_cor_inv <- tryCatch(
    chol2inv(chol(X_reg)),
    error = function(e) solve(X_reg)
  )
  D <- diag(X)
  X_inv <- diag(1 / sqrt(D)) %*% X_cor_inv %*% diag(1 / sqrt(D))

  list(
    inv = X_inv,
    lambda = reg_scale,
    condition = kappa(X_reg, exact = FALSE)
  )
}

.smixcan_col_sds <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colSds(x))
  }
  sqrt(colSums((x - matrix(colMeans(x), nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2) /
    (nrow(x) - 1))
}

.smixcan_col_vars <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::colVars(x))
  }
  colSums((x - matrix(colMeans(x), nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2) /
    (nrow(x) - 1)
}

.estimate_reg_scale_from_cov <- function(X, target_condition = 100,
                                         min_scale = 0.001,
                                         max_scale = 0.1) {
  X_cor <- stats::cov2cor(X)
  eig <- eigen(X_cor, symmetric = TRUE, only.values = TRUE)$values
  lambda_max <- max(eig, na.rm = TRUE)
  lambda_min <- min(eig, na.rm = TRUE)
  if (!is.finite(lambda_max) || !is.finite(lambda_min)) {
    return(max_scale)
  }
  if (lambda_min > 0 && lambda_max / lambda_min <= target_condition) {
    return(min_scale)
  }

  # For X_cor + sI, condition = (lambda_max + s) / (lambda_min + s).
  # Solve for the smallest s that makes condition <= target_condition.
  estimated <- (lambda_max - target_condition * lambda_min) / (target_condition - 1)
  estimated <- max(min_scale, estimated, na.rm = TRUE)
  min(estimated, max_scale)
}

.joint_p_from_inverse <- function(S, v, drop_intercept = FALSE) {
  Z_full <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)
  Z_join <- as.numeric(Z_full)
  if (drop_intercept) {
    Z_join <- Z_join[-1]
  }
  p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)
  list(
    Z_join = Z_join,
    p_join_vec = p_join_vec,
    p_join = safe_ACAT(p_join_vec)
  )
}

#' Run the K-cell-type S-MiXcan association test
#'
#' This function combines cell-type-specific prediction weights, GWAS summary
#' statistics, and a reference genotype panel to compute per-cell-type joint
#' association statistics. The implementation avoids constructing the full
#' SNP-by-SNP LD covariance matrix, reuses predicted expression, and supports
#' fixed or estimated regularization.
#'
#' @name SMiXcan_assoc_test_K
#' @title S-MiXcan Association Test in K cell types
#'
#' @param W Numeric p-by-K matrix of cell-type-specific SNP weights, where
#'   p is the number of SNPs and K is the number of cell types.
#' @param gwas_results List containing GWAS summary statistics with elements
#'   \code{Beta} and \code{se_Beta}.
#' @param x_g Reference-panel genotype matrix with individuals in rows and
#'   SNPs in columns.
#' @param n0 Number of controls.
#' @param n1 Number of cases.
#' @param family Either \code{"binomial"} or \code{"gaussian"}.
#' @param regularization Either \code{"fixed"} or \code{"estimate"}.
#' @param reg_scale Fixed diagonal regularization scale when
#'   \code{regularization = "fixed"}.
#' @param condition_max Target maximum condition number when
#'   \code{regularization = "estimate"}.
#' @param reg_min_scale Minimum estimated regularization scale.
#' @param reg_max_scale Maximum estimated regularization scale.
#' @param weight_eps SNPs with all absolute weights <= \code{weight_eps} are
#'   removed before testing.
#'
#' @return A list containing:
#' \describe{
#'   \item{Z_join}{Vector of per-cell-type joint Z-scores.}
#'   \item{p_join_vec}{Vector of per-cell-type joint p-values.}
#'   \item{p_join}{Combined ACAT p-value across the K cell types.}
#'   \item{Z_sep}{Vector of separate component Z-scores.}
#'   \item{p_sep}{Vector of separate component p-values.}
#'   \item{reg_scale_selected}{Regularization scale used.}
#'   \item{reg_condition}{Condition number after regularization.}
#'   \item{mode}{Either \code{"joint"} or \code{"separate"}.}
#' }
#'
#' @export
SMiXcan_assoc_test_K <- function(W,
                                 gwas_results,
                                 x_g,
                                 n0,
                                 n1,
                                 family = c("binomial", "gaussian"),
                                 regularization = c("fixed", "estimate"),
                                 reg_scale = 0.1,
                                 condition_max = 100,
                                 reg_min_scale = 0.001,
                                 reg_max_scale = 0.1,
                                 weight_eps = 1e-12) {
  family <- match.arg(family)
  regularization <- match.arg(regularization)

  Beta <- as.numeric(gwas_results$Beta)
  se_Beta <- as.numeric(gwas_results$se_Beta)

  if (!is.matrix(W)) {
    stop("W must be a numeric matrix of dimension p * K.")
  }

  p <- nrow(W)
  K <- ncol(W)

  if (ncol(x_g) != p) {
    stop("Number of SNPs (columns) in x_g must match nrow(W).")
  }
  if (length(Beta) != p || length(se_Beta) != p) {
    stop("Length of Beta and se_Beta must match nrow(W).")
  }
  if (any(se_Beta <= 0, na.rm = TRUE)) {
    stop("se_Beta must be positive.")
  }

  sig_l <- .smixcan_col_sds(x_g)
  Yhat <- x_g %*% W
  sig2_g <- .smixcan_col_vars(Yhat)
  num <- colSums(W * (Beta * sig_l / se_Beta))

  Z_sep <- rep(NA_real_, K)
  ok <- sig2_g > 0 & is.finite(sig2_g)
  Z_sep[ok] <- num[ok] / sqrt(sig2_g[ok])
  p_sep <- ifelse(
    is.na(Z_sep),
    NA_real_,
    2 * stats::pnorm(abs(Z_sep), lower.tail = FALSE)
  )

  Z_join <- Z_sep
  p_join_vec <- p_sep
  mode <- "separate"
  reg_scale_selected <- NA_real_
  reg_condition <- NA_real_

  choose_reg_scale <- function(YtY) {
    if (regularization == "estimate") {
      .estimate_reg_scale_from_cov(
        YtY,
        target_condition = condition_max,
        min_scale = reg_min_scale,
        max_scale = reg_max_scale
      )
    } else {
      reg_scale
    }
  }

  if (family == "binomial") {
    if (is.null(n0) || is.null(n1) || is.na(n0) || is.na(n1) || n0 <= 0 || n1 <= 0) {
      stop("n0 and n1 must be positive for binomial family.")
    }

    Z0 <- log(n1 / n0) / sqrt(1 / n1 + 1 / n0)
    if (all(ok)) {
      Y_scaled <- scale(Yhat)
      Y <- cbind(1, Y_scaled)
      colnames(Y) <- c("Y0", paste0("Y", seq_len(K)))
      YtY <- crossprod(Y)
      Omega <- diag(YtY)
      corY <- stats::cov2cor(YtY)
    } else {
      corY <- matrix(NA_real_, nrow = K + 1L, ncol = K + 1L)
    }

    if (any(is.na(corY))) {
      mode <- "invalid_variance_separate"
    } else if (any(abs(corY[upper.tri(corY)]) > 0.999999)) {
      mode <- "collinear_separate"
    } else {
      reg_scale_selected <- choose_reg_scale(YtY)
      inv <- regularized_inverse_cov(YtY, reg_scale = reg_scale_selected)
      v <- c(sqrt(Omega[1]) * Z0, sqrt(Omega[-1]) * Z_sep)
      fit <- .joint_p_from_inverse(inv$inv, v, drop_intercept = TRUE)
      Z_join <- fit$Z_join
      p_join_vec <- fit$p_join_vec
      reg_condition <- inv$condition
      mode <- "joint"
    }
  } else if (family == "gaussian") {
    if (all(ok)) {
      Y_scaled <- scale(Yhat)
      colnames(Y_scaled) <- paste0("Y", seq_len(K))
      YtY <- crossprod(Y_scaled)
      Omega <- diag(YtY)
      corY <- stats::cov2cor(YtY)
    } else {
      corY <- matrix(NA_real_, nrow = K, ncol = K)
    }

    if (any(is.na(corY))) {
      mode <- "invalid_variance_separate"
    } else if (any(abs(corY[upper.tri(corY)]) > 0.999999)) {
      mode <- "collinear_separate"
    } else {
      reg_scale_selected <- choose_reg_scale(YtY)
      inv <- regularized_inverse_cov(YtY, reg_scale = reg_scale_selected)
      v <- sqrt(Omega) * Z_sep
      fit <- .joint_p_from_inverse(inv$inv, v, drop_intercept = FALSE)
      Z_join <- fit$Z_join
      p_join_vec <- fit$p_join_vec
      reg_condition <- inv$condition
      mode <- "joint"
    }
  } else {
    stop("family must be 'binomial' or 'gaussian'")
  }

  p_join <- safe_ACAT(p_join_vec)

  list(
    Z_join = Z_join,
    p_join_vec = p_join_vec,
    p_join = p_join,
    Z_sep = Z_sep,
    p_sep = p_sep,
    reg_scale_selected = reg_scale_selected,
    reg_condition = reg_condition,
    mode = mode
  )
}

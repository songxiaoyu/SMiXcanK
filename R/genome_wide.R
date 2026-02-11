#' Run Single-SNP GWAS (Gaussian or Binomial)
#'
#' Performs univariate GWAS for each SNP in matrix \code{X}.
#' Supports quantitative traits (Gaussian) and case-control traits (Binomial).
#'
#' For Gaussian traits, a fast residualization-based approach is used.
#' For Binomial traits, either per-SNP GLM or a score test can be used.
#'
#' @param X Numeric matrix of genotypes (n samples × p SNPs).
#' @param D Numeric phenotype vector of length n.
#' @param family0 Model family: \code{gaussian()} or \code{binomial()}.
#' @param covar Optional matrix/data.frame of covariates (n × q).
#' @param method_binomial Character. Either \code{"score"} (fast) or \code{"glm"}.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{Beta}{Estimated SNP effect (NA for score test).}
#'   \item{se_Beta}{Standard error (NA for score test).}
#'   \item{Z}{Z-score.}
#'   \item{P}{Two-sided p-value.}
#' }
#'
#' @examples
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' D <- rnorm(n)
#'
#' res <- run_gwas(X, D, gaussian())
#'
#' @export
run_gwas <- function(
    X,
    D,
    family0 = stats::gaussian(),
    covar = NULL,
    method_binomial = c("score", "glm")
) {

  X <- as.matrix(X)
  D <- as.numeric(D)
  n <- nrow(X)
  p <- ncol(X)

  if (length(D) != n) {
    stop("Length of D must equal nrow(X).")
  }

  if (!is.null(covar)) {
    Z <- as.matrix(covar)
    if (nrow(Z) != n) {
      stop("Covariates must have same number of rows as X.")
    }
    Z <- cbind(Intercept = 1, Z)
  } else {
    Z <- matrix(1, n, 1)
    colnames(Z) <- "Intercept"
  }

  is_gaussian <- identical(family0$family, "gaussian")
  is_binomial <- identical(family0$family, "binomial")

  if (!is_gaussian && !is_binomial) {
    stop("Only gaussian() and binomial() are supported.")
  }

  # =============================
  # Gaussian Trait (Fast)
  # =============================
  if (is_gaussian) {

    fit_y <- stats::lm.fit(Z, D)
    ry <- fit_y$residuals

    ZtZ_inv <- solve(crossprod(Z))
    Px <- Z %*% (ZtZ_inv %*% crossprod(Z, X))
    Xr <- X - Px

    XtX <- colSums(Xr * Xr)
    XtY <- as.numeric(crossprod(Xr, ry))

    beta <- XtY / XtX
    df <- n - ncol(Z) - 1
    sigma2 <- sum(ry^2) / max(df, 1)
    se <- sqrt(sigma2 / XtX)

    z <- beta / se
    pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

    return(data.frame(
      Beta = beta,
      se_Beta = se,
      Z = z,
      P = pval
    ))
  }

  # =============================
  # Binomial Trait
  # =============================
  method_binomial <- match.arg(method_binomial)

  if (method_binomial == "glm") {

    results <- matrix(NA_real_, p, 4)

    for (j in seq_len(p)) {

      xj <- X[, j]
      dat <- data.frame(D = D, x = xj)

      if (ncol(Z) > 1) {
        dat <- cbind(dat, Z[, -1, drop = FALSE])
      }

      fit <- tryCatch(
        stats::glm(D ~ ., data = dat, family = family0),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        cs <- summary(fit)$coefficients
        if ("x" %in% rownames(cs)) {
          b <- cs["x", 1]
          se <- cs["x", 2]
          z <- b / se
          pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
          results[j, ] <- c(b, se, z, pval)
        }
      }
    }

    colnames(results) <- c("Beta", "se_Beta", "Z", "P")
    return(as.data.frame(results))
  }

  # =============================
  # Binomial Score Test (Fast)
  # =============================
  null_fit <- stats::glm.fit(x = Z, y = D, family = family0)
  mu <- null_fit$fitted.values
  W <- mu * (1 - mu)
  r <- D - mu

  ZtWZ_inv <- solve(crossprod(Z, Z * W))

  Zmat <- Z
  results <- matrix(NA_real_, p, 2)

  for (j in seq_len(p)) {

    x <- X[, j]
    ZtWx <- crossprod(Zmat, W * x)
    adj <- Zmat %*% (ZtWZ_inv %*% ZtWx)
    xt <- x - adj

    U <- sum(xt * r)
    V <- sum(W * xt * xt)

    if (is.finite(V) && V > 0) {
      z <- U / sqrt(V)
      pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      results[j, ] <- c(z, pval)
    }
  }

  colnames(results) <- c("Z", "P")

  return(data.frame(
    Beta = NA_real_,
    se_Beta = NA_real_,
    Z = results[, 1],
    P = results[, 2]
  ))
}

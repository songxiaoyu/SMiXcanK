#' Safely combine p-values with ACAT
#'
#' This helper wraps \code{ACAT::ACAT()} and returns \code{NA} instead of
#' stopping if ACAT fails because of numerical or input issues.
#'
#' @param p_values Numeric vector of p-values to combine.
#'
#' @return A single combined p-value, or \code{NA} if ACAT fails.
#'
#' @importFrom ACAT ACAT
#' @export
safe_ACAT <- function(p_values) {
  p_values <- as.numeric(p_values)
  p_values <- p_values[is.finite(p_values) & p_values >= 0 & p_values <= 1]
  if (length(p_values) == 0L) {
    return(NA_real_)
  }
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)
  })
}


#' Cell-type-aware association test for two predicted expression traits
#'
#' This function performs a two-cell-type association analysis using predicted
#' expression values from two cell types. When the two predictors are nearly
#' collinear, it falls back to separate univariate models; otherwise it
#' applies the shrinkage-based joint test used by MiXcan.
#'
#' @param outcome Numeric or factor phenotype vector.
#' @param cell1 Numeric vector of predicted expression for cell type 1.
#' @param cell2 Numeric vector of predicted expression for cell type 2.
#' @param family Either \code{"binomial"} or \code{"gaussian"}.
#' @param rho_thr Correlation threshold above which the function uses the
#'   separate-model fallback.
#'
#' @return A list containing cell-type-specific effect estimates, standard
#'   errors, p-values, a combined ACAT p-value, and the fitting mode used.
#' @importFrom stats model.matrix glm binomial na.omit
#' @importFrom dplyr select bind_cols mutate
#' @export
MiXcan_assoc_test <- function(outcome, cell1, cell2, family = 'binomial',
                                   rho_thr = 0.999999) {
  # Coerce inputs to numeric vectors.
  y  <- as.numeric(outcome)
  y1 <- as.numeric(cell1)
  y2 <- as.numeric(cell2)

  df <- data.frame(outcome = y, Y1 = y1, Y2 = y2)
  df <- stats::na.omit(df)
  if (nrow(df) == 0L) {
    return(list(cell1_est=NA, cell1_se=NA, cell1_p=NA,
                cell2_est=NA, cell2_se=NA, cell2_p=NA,
                p_combined=NA, mode="empty"))
  }

  # Standardize the predicted-expression inputs.
  Ys <- scale(as.matrix(df[, c("Y1","Y2")]), center = TRUE, scale = TRUE)
  colnames(Ys) <- c("Y1","Y2")

  # Null intercept-only model used to form the baseline z-score.
  fit0 <- glm(df$outcome ~ 1, family = family)
  coef0 <- summary(fit0)$coefficients
  if (family == "gaussian") {
    Z0 <- coef0["(Intercept)", "t value"]
  } else {
    Z0 <- coef0["(Intercept)", "z value"]
  }

  # Fit separate one-predictor models for the fallback path.
  f1 <- glm(df$outcome ~ Ys[, "Y1"], family = family)
  f2 <- glm(df$outcome ~ Ys[, "Y2"], family = family)
  s1 <- summary(f1)$coefficients
  s2 <- summary(f2)$coefficients

  # Extract estimates, standard errors, z-scores, and p-values safely.
  get_stats <- function(s, idx = 2L) {
    if (nrow(s) >= 2L && all(c("Estimate","Std. Error") %in% colnames(s))) {
      est <- s[idx, "Estimate"]; se <- s[idx, "Std. Error"]
      z   <- est / se
      p   <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      c(est = est, se = se, z = z, p = p)
    } else c(est = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_)
  }
  g1 <- get_stats(s1); g2 <- get_stats(s2)

  # Check whether the two predictors are nearly collinear.
  rho <- suppressWarnings(stats::cor(Ys[,1], Ys[,2]))
  if (!is.finite(rho)) rho <- 0

  # Fallback to separate models if the design is nearly singular.
  if (abs(rho) > rho_thr) {
    p_comb <- safe_ACAT(c(g1["p"], g2["p"]))
    return(list(
      cell1_est = unname(g1["est"]), cell1_se = unname(g1["se"]), cell1_p = unname(g1["p"]),
      cell2_est = unname(g2["est"]), cell2_se = unname(g2["se"]), cell2_p = unname(g2["p"]),
      p_combined = p_comb, mode = "collinear_separate"
    ))
  }
  Y  <- cbind(1, Ys)                 # n×3 : intercept, Y1, Y2
  YtY <- crossprod(Y)                # 3×3
  Omega <- diag(YtY)


  S <- regularized_inverse_cov(YtY)$inv
  v <- c(sqrt(Omega[1]) * Z0, sqrt(Omega[2]) * g1["z"], sqrt(Omega[3]) * g2["z"])
  Z_join <- diag(1 / sqrt(diag(S))) %*% S %*% matrix(v, ncol = 1)

  Z1_join <- as.numeric(Z_join[2])
  Z2_join <- as.numeric(Z_join[3])
  p1 <- 2 * stats::pnorm(abs(Z1_join), lower.tail = FALSE)
  p2 <- 2 * stats::pnorm(abs(Z2_join), lower.tail = FALSE)
  p_comb <- safe_ACAT(c(p1, p2))

  # Return estimates and SEs from the joint two-predictor fit.
  fit_biv <- glm(df$outcome ~ Ys, family = family)
  sb <- summary(fit_biv)$coefficients
  est1 <- if ("YsY1" %in% rownames(sb)) sb["YsY1","Estimate"] else NA_real_
  se1  <- if ("YsY1" %in% rownames(sb)) sb["YsY1","Std. Error"] else NA_real_
  est2 <- if ("YsY2" %in% rownames(sb)) sb["YsY2","Estimate"] else NA_real_
  se2  <- if ("YsY2" %in% rownames(sb)) sb["YsY2","Std. Error"] else NA_real_

  list(
    cell1_est = est1, cell1_se = se1, cell1_p = p1,
    cell2_est = est2, cell2_se = se2, cell2_p = p2,
    p_combined = p_comb, mode = "joint"
  )
}

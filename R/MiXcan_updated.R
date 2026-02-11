#' Safe Wrapper for ACAT P-value Combination
#'
#' This function computes a combined p-value using the ACAT method while handling possible numerical or input errors.
#' If an error occurs during computation, the function returns \code{NA} and prints an informative message.
#'
#' @param p_values A numeric vector of p-values to be combined using the ACAT method.
#'
#' @return A single numeric value: the combined p-value (or \code{NA} if ACAT fails).
#'
#' @importFrom ACAT ACAT
#' @export
safe_ACAT <- function(p_values) {
  tryCatch({
    ACAT::ACAT(p_values)
  }, error = function(e) {
    message("ACAT Error occurred: ", e$message)
    return(NA)
  })
}


#' @title Cell-Type-Aware Association with Shrinkage
#'
#' @description
#' This function performs cell-type-aware association analysis with penalization. It estimates the individual and combined p-values for
#' predicted gene expressions (e.g., from cell type 1 and cell type 2).
#' @param outcome A binary phenotype vector (0/1 or factor).
#' @param cell1 y predicted for cell type1.
#' @param cell2 y predicted for cell type2.
#' @param family 'binomial' or 'gaussian'
#' @param rho_thr threshold where choose sep model
#' @return A data frame containing effect size estimates, standard errors, and p-values for both cell types,
#'         and a combined p-value from ACAT.
#' @importFrom stats model.matrix glm binomial na.omit
#' @importFrom dplyr select bind_cols mutate
#' @export
MiXcan_assoc_test <- function(outcome, cell1, cell2, family = 'binomial',
                                   rho_thr = 0.999999) {
  # coerce to numeric vectors (IMPORTANT: outcome as vector, not matrix)
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

  # standardize predictors
  Ys <- scale(as.matrix(df[, c("Y1","Y2")]), center = TRUE, scale = TRUE)
  colnames(Ys) <- c("Y1","Y2")

  # null/intercept-only logistic => Z0
  fit0 <- glm(df$outcome ~ 1, family = family)
  if (family == "gaussian") {
    Z0 <- coef0["(Intercept)", "t value"]
  } else {
    Z0 <- coef0["(Intercept)", "z value"]
  }

  # separate univariate logistic z/p
  f1 <- glm(df$outcome ~ Ys[, "Y1"], family = family)
  f2 <- glm(df$outcome ~ Ys[, "Y2"], family = family)
  s1 <- summary(f1)$coefficients
  s2 <- summary(f2)$coefficients

  # robust extraction of z/est/se
  get_stats <- function(s, idx = 2L) {
    if (nrow(s) >= 2L && all(c("Estimate","Std. Error") %in% colnames(s))) {
      est <- s[idx, "Estimate"]; se <- s[idx, "Std. Error"]
      z   <- est / se
      p   <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      c(est = est, se = se, z = z, p = p)
    } else c(est = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_)
  }
  g1 <- get_stats(s1); g2 <- get_stats(s2)

  # correlation check
  rho <- suppressWarnings(stats::cor(Ys[,1], Ys[,2]))
  if (!is.finite(rho)) rho <- 0

  # (1) Fallback if near-singular
  if (abs(rho) > rho_thr) {
    p_comb <- safe_ACAT(c(g1["p"], g2["p"]))
    return(list(
      cell1_est = unname(g1["est"]), cell1_se = unname(g1["se"]), cell1_p = unname(g1["p"]),
      cell2_est = unname(g2["est"]), cell2_se = unname(g2["se"]), cell2_p = unname(g2["p"]),
      p_combined = p_comb, mode = "separate"
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

  # also return estimates/SEs from the bivariate fit (optional)
  fit_biv <- glm(df$outcome ~ Ys, family = family)
  sb <- summary(fit_biv)$coefficients
  est1 <- if ("YsY1" %in% rownames(sb)) sb["YsY1","Estimate"] else NA_real_
  se1  <- if ("YsY1" %in% rownames(sb)) sb["YsY1","Std. Error"] else NA_real_
  est2 <- if ("YsY2" %in% rownames(sb)) sb["YsY2","Estimate"] else NA_real_
  se2  <- if ("YsY2" %in% rownames(sb)) sb["YsY2","Std. Error"] else NA_real_

  list(
    cell1_est = est1, cell1_se = se1, cell1_p = p1,
    cell2_est = est2, cell2_se = se2, cell2_p = p2,
    p_combined = p_comb
  )
}


#' Estimate cell-type-level prediction weights for a gene (symmetric K-cell MiXcan)
#'
#' This function generalizes the original MiXcan two–cell-type model to
#' K cell types using a symmetric mean + contrast parameterization.
#' Gene expression is modeled as a linear combination of:
#'   (i) shared SNP effects across all cell types,
#'   (ii) contrast-coded cell-type effects, and
#'   (iii) SNP-by-contrast interactions.
#'
#' Elastic-net regression is used with no penalty on the contrast-only
#' terms corresponding to differences in baseline expression across cell
#' types, and equal penalization on shared SNP effects and contrast SNP
#' effects, ensuring symmetric treatment of all cell types.
#'
#' @param y Numeric vector or N × 1 matrix.
#'   Pre-cleaned expression levels for a single gene in N samples.
#' @param x Numeric N × P matrix.
#'   Genotype matrix for cis-SNP predictors of the gene.
#' @param cov Optional N × Q matrix of covariates (e.g., age, PCs).
#'   If \code{NULL}, no covariates are included. Default: \code{NULL}.
#' @param pi_k Numeric N × K matrix of estimated cell-type fractions for
#'   K cell types. Each row corresponds to a sample, each column to a
#'   cell type. Rows are typically approximately summing to 1.
#' @param xNameMatrix Optional SNP annotation (e.g., rsID, position, etc.).
#'   If \code{NULL}, generic names \code{"SNP1", ..., "SNPp"} are used.
#' @param yName Optional gene annotation (e.g., gene ID, symbol).
#' @param foldid Optional length-N vector of CV fold IDs. If \code{NULL},
#'   10-fold CV with a random split is used.
#' @param alpha Elastic-net mixing parameter passed to \code{glmnet}.
#'   Default: \code{0.5}.
#'
#' @return A list with components:
#' \describe{
#'   \item{type}{Character string: \code{"CellTypeSpecific"}, \code{"NonSpecific"},
#'     or \code{"NoPredictor"}, based on whether contrast SNP effects are present.}
#'   \item{beta.SNP.by.cell}{List of length K; each element is a data.frame
#'     containing SNP annotation and cell-type-specific SNP weights.}
#'   \item{beta.all.models}{Matrix of regression coefficients for the tissue-
#'     level model and K cell-type models (intercepts + SNP + covariates).}
#'   \item{W}{Numeric matrix (SNPs × K) of cell-type-specific SNP weights.
#'     This is the matrix used directly in \code{SMiXcan_assoc_test_K()}.}
#'   \item{glmnet.cell}{Fitted glmnet object for the K-cell model.}
#'   \item{glmnet.tissue}{Fitted glmnet object for the tissue-level model
#'     (ignoring cell-type composition).}
#'   \item{yName}{Returned gene annotation.}
#'   \item{xNameMatrix}{Returned SNP annotation.}
#' }
#'
#' @details
#' The model uses a contrast matrix \eqn{C} (by default Helmert contrasts)
#' of size K × (K-1). For sample i with cell-type fractions \eqn{\pi_{ik}},
#' contrast-coded compositions are defined as
#' \deqn{c_i = \pi_i^\top C}
#' and SNP-by-contrast interactions are constructed as
#' \deqn{(c_{im} x_{il}), \quad m=1,\dots,K-1; \; l=1,\dots,P.}
#'
#' The linear predictor is
#' \deqn{
#'   y_i = \beta_0 + \sum_{m=1}^{K-1} c_{im} \alpha_m
#'         + x_i^\top \bar{b}
#'         + \sum_{m=1}^{K-1} c_{im} x_i^\top d_m
#'         + z_i^\top \gamma + \varepsilon_i,
#' }
#' where \eqn{\bar{b}} is the shared SNP effect, \eqn{d_m} are contrast
#' effects, and \eqn{z_i} are optional covariates.
#' Cell-type-specific SNP effects are reconstructed as
#' \deqn{
#'   b_k = \bar{b} + \sum_{m=1}^{K-1} C_{km} d_m,
#' }
#' ensuring symmetric treatment of all cell types. For K=2, this reduces
#' to the original MiXcan parameterization up to a constant scaling of
#' the contrast.
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @export
MiXcan_train_K <- function(y, x, cov = NULL, pi_k,
                                     xNameMatrix = NULL, yName = NULL,
                                     foldid = NULL, alpha = 0.5) {

  ## ---- Basic checks and setup ----
  y    <- as.matrix(y)
  x    <- as.matrix(x)
  pi_k <- as.matrix(pi_k)

  n <- nrow(x)
  p <- ncol(x)
  K <- ncol(pi_k)

  if (nrow(pi_k) != n) {
    stop("pi_k must have the same number of rows as x (samples).")
  }

  # SNP names / annotation
  if (is.null(xNameMatrix)) {
    if (!is.null(colnames(x))) {
      xNameMatrix <- data.frame(SNP = colnames(x), stringsAsFactors = FALSE)
    } else {
      xNameMatrix <- data.frame(SNP = paste0("SNP", seq_len(p)),
                                stringsAsFactors = FALSE)
    }
  }

  if (is.null(foldid)) {
    set.seed(1L)
    foldid <- sample(1:10, n, replace = TRUE)
  }

  # center y and x for numerical stability; pi_k left on original scale
  y <- scale(y, center = TRUE, scale = FALSE)
  x <- scale(x, center = TRUE, scale = FALSE)

  ## ---- Tissue-level model (PrediXcan-like) ----
  if (is.null(cov)) {
    xcov <- x
    pcov <- 0L
  } else {
    cov  <- as.matrix(cov)
    cov  <- scale(cov, center = TRUE, scale = FALSE)
    pcov <- ncol(cov)
    xcov <- cbind(x, cov)
  }

  ft00 <- glmnet::cv.glmnet(
    x = xcov, y = y,
    family   = "gaussian",
    foldid   = foldid,
    alpha    = alpha
  )
  ft0 <- glmnet::glmnet(
    x = xcov, y = y,
    family   = "gaussian",
    lambda   = ft00$lambda.1se,
    alpha    = alpha
  )
  est.tissue <- c(ft0$a0, as.numeric(ft0$beta))  # intercept + coefficients

  ## ---- Symmetric K-cell model: mean + contrasts ----

  # Contrast matrix C: K x (K-1), e.g. Helmert contrasts
  C <- stats::contr.helmert(K)   # each row corresponds to 1 cell type
  # Contrast-coded compositions: N x (K-1)
  c_mat <- pi_k %*% C

  # Build interaction blocks: for each contrast m, c_mat[,m] * x (elementwise)
  Z_list <- vector("list", length = K - 1L)
  for (m in seq_len(K - 1L)) {
    Z_list[[m]] <- c_mat[, m] * x   # N x P
  }
  Z_block <- do.call(cbind, Z_list)  # N x [P*(K-1)]

  if (is.null(cov)) {
    XX <- cbind(c_mat, x, Z_block)  # N x [ (K-1) + P + P*(K-1) ]
  } else {
    XX <- cbind(c_mat, x, Z_block, cov)
  }

  n_c  <- K - 1L
  n_x  <- p
  n_Z  <- p * (K - 1L)

  # Penalization:
  #  - no penalty on c_mat (contrast-only terms for intercept differences)
  #  - penalty on shared SNP effects, contrast SNP effects, and covariates
  penalty.factor <- c(
    rep(0, n_c),                 # c_mat (unpenalized)
    rep(1, n_x + n_Z + pcov)     # x, Z_block, cov
  )

  ft11 <- glmnet::cv.glmnet(
    x = XX, y = y,
    family         = "gaussian",
    foldid         = foldid,
    alpha          = alpha,
    penalty.factor = penalty.factor
  )
  ft <- glmnet::glmnet(
    x = XX, y = y,
    family         = "gaussian",
    lambda         = ft11$lambda.1se,
    alpha          = alpha,
    penalty.factor = penalty.factor
  )

  est <- c(ft$a0, as.numeric(ft$beta))  # [intercept, all coefficients]

  ## ---- Decode parameters: intercepts and SNP weights per cell type ----

  a0 <- ft$a0   # scalar intercept

  # indices in est (excluding intercept a0):
  # 1..n_c         : alpha_m (effects of c_mat on intercept)
  # (n_c+1)..(n_c+p) : shared SNP effects (bar{b})
  # next n_Z        : contrast SNP effects stacked
  # last pcov       : covariate effects (shared across cell types)

  idx_c    <- 1:n_c
  idx_barb <- (n_c + 1):(n_c + p)
  idx_Z    <- (n_c + p + 1):(n_c + p + n_Z)
  idx_cov  <- if (pcov > 0L) (length(est) - pcov + 1):length(est) else integer(0)

  alpha_vec <- est[idx_c]            # length K-1
  b_bar     <- est[idx_barb]         # length P

  # contrast SNP effects d_m: P x (K-1)
  d_mat <- matrix(est[idx_Z], nrow = p, ncol = K - 1L, byrow = FALSE)

  # reconstruct cell-type-specific intercepts:
  # a_k = a0 + C[k,] %*% alpha_vec
  a_cell <- as.numeric(a0 + C %*% alpha_vec)  # length K

  # reconstruct cell-type-specific SNP weights:
  # b_k = b_bar + d_mat %*% C[k,]
  B_mat <- matrix(NA_real_, nrow = p, ncol = K)
  for (k in seq_len(K)) {
    B_mat[, k] <- b_bar + d_mat %*% C[k, ]
  }
  colnames(B_mat) <- paste0("Cell", seq_len(K))
  rownames(B_mat) <- xNameMatrix[[1]]

  ## ---- Determine model type ----
  # If all contrast SNP effects are exactly zero => NonSpecific / NoPredictor
  if (suppressWarnings(all(d_mat == 0))) {
    if (suppressWarnings(all(b_bar == 0))) {
      Type <- "NoPredictor"
    } else {
      Type <- "NonSpecific"
    }
  } else {
    Type <- "CellTypeSpecific"
  }

  ## ---- Assemble outputs ----

  # 1) per-cell-type SNP weights as data.frames
  beta.SNP.by.cell <- vector("list", length = K)
  for (k in seq_len(K)) {
    beta.SNP.by.cell[[k]] <- data.frame(
      weight = B_mat[, k],
      stringsAsFactors = FALSE
    )
    rownames(beta.SNP.by.cell[[k]]) <- rownames(B_mat)
  }
  names(beta.SNP.by.cell) <- paste0("Cell", seq_len(K))

  # 2) beta.all.models:
  #    combine tissue model + K cell-type models (intercepts + SNP [+ cov])
  n_rows <- 1 + p + pcov
  beta_cell_mat <- matrix(NA_real_, nrow = n_rows, ncol = K)
  rownames(beta_cell_mat) <- c(
    "Intercept",
    xNameMatrix[[1]],
    if (pcov > 0) paste0("cov", seq_len(pcov)) else NULL
  )
  colnames(beta_cell_mat) <- paste0("Cell", seq_len(K))

  for (k in seq_len(K)) {
    beta_cell_mat[1, k] <- a_cell[k]          # intercept
    beta_cell_mat[1 + (1:p), k] <- B_mat[, k] # SNPs
    if (pcov > 0L) {
      beta_cell_mat[1 + p + (1:pcov), k] <-
        if (length(idx_cov) > 0L) est[idx_cov] else 0
    }
  }

  # tissue model: est.tissue already contains [intercept, SNPs, covs]
  beta_all_models <- cbind(
    Tissue = est.tissue,
    beta_cell_mat
  )

  list(
    type             = Type,
    beta.SNP.by.cell = beta.SNP.by.cell,
    beta.all.models  = beta_all_models,
    W                = B_mat,          # ← ADD THIS
    glmnet.cell      = ft,
    glmnet.tissue    = ft0,
    yName            = yName,
    xNameMatrix      = xNameMatrix
  )

}

# --- Stage 2: run_s_mixcan_analysis.R ---

library(data.table)
library(lme4)
library(glmnet)
library(doParallel)
library(doRNG)
library(ACAT)
library(SMiXcan)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)
library(SMiXcanK)

# STEP1: load data--------
setwd('/Users/zhusinan/Downloads/S-MiXcan_code_folder/2pi')

n1 <- 133384
# n1<- 106278
n0 <- 113789 + 18908
# n0 <- 91477

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
    p_join    = p_join
  )
}


# Set chromosome
for(chr in 2:22){
  print('CHR')
  print(chr)
  # Load pre-merged data
  mw_gwas_input_path <- file.path(sprintf("bcac2020_input/chr%d_mw_gwas_input_bcac2020_pi2.rds", chr))
  mw_gwas_input <- readRDS(mw_gwas_input_path)

  # Load LD matrix, ref genome, SNP info
  LD_input_dir <- 'bcac2020_filtered_id/'
  LD_snp_path <- file.path(LD_input_dir,sprintf("filtered_chr%d_hg38_pi2.bim", chr))
  X_ref_path <- file.path(LD_input_dir,sprintf("filtered_chr%d_hg38_012_pi2.raw", chr))

  ld_snp <- fread(LD_snp_path, header = FALSE)
  ref_snp_id <- ld_snp$V2

  # Read genotype matrix and set proper colnames
  X_ref <- as.matrix(fread(X_ref_path)[, 7:ncol(fread(X_ref_path))])
  colnames(X_ref) <- sub("_.*", "", colnames(X_ref))

  # Filter intersecting SNPs
  nrow(mw_gwas_input)
  gwas_ref_snps <- intersect(mw_gwas_input$varID, ref_snp_id)
  length(gwas_ref_snps)
  filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$varID %in% gwas_ref_snps, ]

  # Prepare filtered gene list
  split_df <- split(filtered_mw_gwas_input, filtered_mw_gwas_input$gene)
  filtered_list <- list()
  for (gene in names(split_df)) {
    gene_df <- split_df[[gene]]
    W1 <- gene_df$weight_cell_1
    W2 <- gene_df$weight_cell_2
    W <- cbind(W1, W2, W3)
    filtered_list[[gene]] <- list(W = W, selected_snp = gene_df)
  }


  # STEP2 Run S-MiXcan----
  G <- length(filtered_list)
  real_result = data.frame(matrix(ncol = 10, nrow = G))
  colnames(real_result) <- c('gene_name','gene_id','chr','type','input_snp_num','Z_1','p_1','Z_2','p_2','p_join')
  real_result$chr <- chr

  for (g in 1:G) {
    gene = names(split_df)[g]
    cat("Processing gene:", gene, "\n")

    W <- filtered_list[[gene]]$W
    selected_snp_id <- filtered_list[[gene]]$selected_snp$varID
    gwas_results <- list(
      Beta = filtered_list[[gene]]$selected_snp$beta.Gwas,
      se_Beta = filtered_list[[gene]]$selected_snp$SE.Gwas
    )
    X_ref_filtered <- X_ref[, selected_snp_id, drop = FALSE]
    S_MiXcan_results <- SMiXcan_assoc_test_K(W, gwas_results, X_ref_filtered, n0=n0, n1=n1, family='binomial')
    real_result[g, c('gene_name','gene_id','chr', 'type', 'input_snp_num')] =c(filtered_list[[gene]]$selected_snp[1,c('gene','varID','CHR','type')], nrow(W))
    real_result[g, c('Z_1','p_1','Z_2','p_2','p_join')] <- c(c(rbind(S_MiXcan_results$Z_join, S_MiXcan_results$p_join_vec)), S_MiXcan_results$p_join)
  }

  result_path <- file.path('bcac2020_result',sprintf("bcac2020_chr%d_result_pi2.csv", chr))
  write.csv(real_result, result_path, row.names = FALSE)
}

combined <- do.call(rbind, lapply(1:22, function(chr) {
  result_path <- file.path("bcac2020_result",
                           sprintf("bcac2020_chr%d_result_pi2.csv", chr))
  read.csv(result_path)
}))
combined_path <- file.path('bcac2020_result','bcac2020_result_pi2.csv')
write.csv(combined, combined_path, row.names = FALSE)
# Extract p-values
pvals <- real_result$p_join
pvals <- combined$p_join
#n=6000
#pvals=matrix(runif(3*n, 0,1), ncol=3) # change to your p-values from cell type 1, 2, 3

pvals = combined[,c('p_1','p_2','p_3')]
library("Primo")
#Primo_pval(pvals=pvals, alt_props=c(0.005, 0.005, 0.005))

library(bacon)
merged_1 <- merged[which(merged$type_ct2 == 'CellTypeSpecific')]
pvals <- merged$p_1_ct2[which(merged_1$p_2_ct2 < 0.01)]

# keep p strictly in (0,1) to avoid Inf
p <- pmin(pmax(pvals, .Machine$double.eps), 1 - .Machine$double.eps)
# two-sided Z (unsigned; direction not needed for bacon)
y <- qnorm(1-p, lower.tail = FALSE)

# run bacon, excluding any non-finite values
bc <- bacon(y, na.exclude = TRUE)

pvals <- pval(bc)
# inspect bias and inflation (λ)
estimates(bc)      # bias and inflation
lambda <- inflation(bc)      # same λ
lambda

pvals <- merged$p_1_ct2
# Expected vs Observed
expected <- -log10(ppoints(length(pvals)))
observed <- -log10(sort(pvals))

# run genome inflation
# Genomic inflation factor
#lambda <- 1.007253

# Save
#pdf("QQplot of .pdf", width=5, height=5)

# QQ plot
plot(expected, observed,
     xlab=expression(Expected~~-log[10](italic(p))),
     ylab=expression(Observed~~-log[10](italic(p))),
     main="QQ Plot of SMiXcan_K in 2 cell types",
     pch=19, cex=0.6, col="blue", las=1)
abline(0, 1, col="red", lwd=2, lty=2)

text(x=min(expected), y=max(observed)*0.9,
     labels=bquote(lambda["GC"] == .(round(lambda,3))),
     adj=0, cex=1)

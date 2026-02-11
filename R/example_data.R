# ==============================================================
# Strong-signal example data for SMiXcanK (2 cell types)
# Goal: MiXcan_train_K() returns "CellTypeSpecific" stably
# ==============================================================

set.seed(2026)

# ----------------------------
# 1) Dimensions
# ----------------------------
n <- 20
p <- 10
K <- 2

# ----------------------------
# 2) Genotype matrix X (N x P)
#    Use higher MAF so SNPs have variance
# ----------------------------
x_example <- matrix(
  rbinom(n * p, 2, 0.5),   # MAF ~ 0.5
  nrow = n,
  ncol = p
)
colnames(x_example) <- paste0("SNP", 1:p)
rownames(x_example) <- paste0("Sample", 1:n)

# ----------------------------
# 3) Cell-type fractions pi_k (N x 2)
#    Make near-pure compositions (extreme) to strengthen contrasts
# ----------------------------
pi_k <- matrix(NA_real_, nrow = n, ncol = K)

# First half: almost all Cell1
pi_k[1:(n/2), 1] <- runif(n/2, 0.97, 0.995)
pi_k[1:(n/2), 2] <- 1 - pi_k[1:(n/2), 1]

# Second half: almost all Cell2
pi_k[(n/2 + 1):n, 2] <- runif(n/2, 0.97, 0.995)
pi_k[(n/2 + 1):n, 1] <- 1 - pi_k[(n/2 + 1):n, 2]

colnames(pi_k) <- c("Cell1", "Cell2")
rownames(pi_k) <- rownames(x_example)

# ----------------------------
# 4) Construct *very strong* cell-type-specific SNP effects
#    Make Cell1 and Cell2 effects opposite for first 6 SNPs
# ----------------------------
b1 <- c( 8, -8,  6, -6,  5, -5, 0, 0, 0, 0)  # Cell1
b2 <- c(-8,  8, -6,  6, -5,  5, 0, 0, 0, 0)  # Cell2 (opposite)

# ----------------------------
# 5) Generate expression y with strong contrast signal
#    y_i = pi1_i * (X b1) + pi2_i * (X b2) + small noise
# ----------------------------
signal1 <- as.numeric(x_example %*% b1)
signal2 <- as.numeric(x_example %*% b2)

y_example <- pi_k[, 1] * signal1 +
  pi_k[, 2] * signal2 +
  rnorm(n, sd = 0.05)   # very small noise

names(y_example) <- rownames(x_example)

# ----------------------------
# 6) Minimal GWAS summary stats (length p)
#    (Just demo format; not meant for real inference)
# ----------------------------
gwas_example <- list(
  Beta    = rnorm(p, 0, 0.05),
  se_Beta = rep(0.05, p)
)
names(gwas_example$Beta)    <- colnames(x_example)
names(gwas_example$se_Beta) <- colnames(x_example)

## -----------------------------
## 7) Example genome-wide results
## -----------------------------

set.seed(123)

G <- 300
gene <- paste0("Gene", seq_len(G))

# 150 specific, 150 nonspecific
type_ct2 <- c(rep("CellTypeSpecific", G/2),
              rep("NonSpecific",     G/2))

# baseline null p-values
p_1_ct2    <- runif(G)
p_2_ct2    <- runif(G)
p_join_ct2 <- runif(G)

## ---- Inject true cell-type-specific signals ----
idx_spec <- which(type_ct2 == "CellTypeSpecific")

# Cell1-only
sig1 <- sample(idx_spec, 20)
p_1_ct2[sig1] <- 10^runif(length(sig1), -8, -5)

# Cell2-only
remaining <- setdiff(idx_spec, sig1)
sig2 <- sample(remaining, 20)
p_2_ct2[sig2] <- 10^runif(length(sig2), -8, -5)

# Shared specific
remaining <- setdiff(idx_spec, c(sig1, sig2))
shared_spec <- sample(remaining, 15)
p_1_ct2[shared_spec] <- 10^runif(length(shared_spec), -8, -5)
p_2_ct2[shared_spec] <- 10^runif(length(shared_spec), -8, -5)

## ---- Inject nonspecific joint signals ----
idx_uns <- which(type_ct2 == "NonSpecific")
sig_joint <- sample(idx_uns, 25)
p_join_ct2[sig_joint] <- 10^runif(length(sig_joint), -10, -6)

merged_example <- data.frame(
  gene        = gene,
  type_ct2    = type_ct2,
  p_1_ct2     = p_1_ct2,
  p_2_ct2     = p_2_ct2,
  p_join_ct2  = p_join_ct2,
  stringsAsFactors = FALSE
)

# ----------------------------
# 8) Save into package data/
# ----------------------------
usethis::use_data(
  x_example, y_example, pi_k, gwas_example, merged_example,
  overwrite = TRUE
)

cat("Saved strong-signal example datasets: x_example, y_example, pi_k, gwas_example, merged_example\n")

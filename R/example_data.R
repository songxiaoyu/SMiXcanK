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

# ----------------------------
# 7) merged_example for PRIMO demo
#    Make sure there are enough non-null-ish signals to avoid Primo failures
# ----------------------------
merged_example <- data.frame(
  gene = paste0("Gene", 1:80),
  type_ct2 = c(rep("CellTypeSpecific", 50), rep("NonSpecific", 30)),
  p_1_ct2 = c(runif(10, 1e-10, 1e-4), runif(40, 0.2, 1), runif(30, 0.2, 1)),
  p_2_ct2 = c(runif(10, 0.2, 1),      runif(10, 1e-10, 1e-4), runif(60, 0.2, 1)),
  p_join_ct2 = c(runif(50, 0.2, 1), runif(5, 1e-12, 1e-6), runif(25, 0.2, 1))
)

# ----------------------------
# 8) Save into package data/
# ----------------------------
usethis::use_data(
  x_example, y_example, pi_k, gwas_example, merged_example,
  overwrite = TRUE
)

cat("Saved strong-signal example datasets: x_example, y_example, pi_k, gwas_example, merged_example\n")

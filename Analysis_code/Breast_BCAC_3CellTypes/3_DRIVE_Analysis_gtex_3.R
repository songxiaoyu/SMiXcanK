# Validate the three-cell-type MiXcan model in the individual-level DRIVE data.
#
# DRIVE genotype data currently cover chromosome 21. This script compares
# individual-level MiXcan p-values with S-MiXcan p-values calculated from a
# single-SNP GWAS in the same DRIVE samples. The three modeled components are:
#   1. Adipocyte/endothelial
#   2. Fibroblast
#   3. Epithelial

suppressPackageStartupMessages({
  library(data.table)
  library(SMiXcan)
})

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
drive_dir <- file.path(paper_dir, "Data", "DRIVE")
results_dir <- file.path(paper_dir, "Results")

weight_path <- file.path(results_dir, "weights_miXcan_full_pi3.csv")
genotype_path <- file.path(drive_dir, "oncoarray_dosages_chr21.txt")
phenotype_path <- file.path(drive_dir, "oncoarray-drive.pheno.csv")
output_path <- file.path(results_dir, "drive_result_pi3.csv")

required_files <- c(weight_path, genotype_path, phenotype_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

# Generalize the regularized individual-level MiXcan association calculation
# to K predicted-expression components. This mirrors MiXcan_assoc_test() while
# returning vectors for all three cell types.
mixcan_assoc_test_k <- function(outcome, predicted_expression,
                                family = "binomial", reg_scale = 0.1,
                                correlation_limit = 0.999999) {
  y <- as.numeric(outcome)
  Y <- as.matrix(predicted_expression)
  keep <- complete.cases(y, Y)
  y <- y[keep]
  Y <- Y[keep, , drop = FALSE]
  K <- ncol(Y)

  empty <- list(
    p_sep = rep(NA_real_, K),
    p_sep_combined = NA_real_,
    p_join_vec = rep(NA_real_, K),
    p_join = NA_real_,
    mode = "empty"
  )
  if (length(y) == 0L || length(unique(y)) < 2L || nrow(Y) <= K + 1L) {
    return(empty)
  }

  Y_scaled <- scale(Y)
  if (any(!is.finite(Y_scaled))) {
    return(empty)
  }
  colnames(Y_scaled) <- paste0("Y", seq_len(K))

  fit0 <- tryCatch(
    stats::glm(y ~ 1, family = family),
    error = function(e) NULL
  )
  if (is.null(fit0)) {
    return(empty)
  }
  coef0 <- summary(fit0)$coefficients
  stat_col <- if (family == "gaussian") "t value" else "z value"
  Z0 <- coef0["(Intercept)", stat_col]

  Z_sep <- vapply(seq_len(K), function(k) {
    fit <- tryCatch(
      stats::glm(y ~ Y_scaled[, k], family = family),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NA_real_)
    }
    cs <- summary(fit)$coefficients
    if (nrow(cs) < 2L || !all(c("Estimate", "Std. Error") %in% colnames(cs))) {
      return(NA_real_)
    }
    cs[2L, "Estimate"] / cs[2L, "Std. Error"]
  }, numeric(1))
  p_sep <- 2 * stats::pnorm(abs(Z_sep), lower.tail = FALSE)

  result <- list(
    p_sep = p_sep,
    p_sep_combined = SMiXcan::safe_ACAT(p_sep),
    p_join_vec = p_sep,
    p_join = SMiXcan::safe_ACAT(p_sep),
    mode = "separate"
  )

  design <- cbind(Intercept = 1, Y_scaled)
  design_crossprod <- crossprod(design)
  design_cor <- suppressWarnings(stats::cov2cor(design_crossprod))
  if (any(!is.finite(design_cor)) ||
      any(abs(design_cor[upper.tri(design_cor)]) > correlation_limit)) {
    return(result)
  }

  inverse <- SMiXcan::regularized_inverse_cov(
    design_crossprod,
    reg_scale = reg_scale
  )$inv
  omega <- diag(design_crossprod)
  score <- c(sqrt(omega[1L]) * Z0, sqrt(omega[-1L]) * Z_sep)
  Z_join <- as.numeric(
    diag(1 / sqrt(diag(inverse))) %*% inverse %*% matrix(score, ncol = 1)
  )[-1L]
  p_join_vec <- 2 * stats::pnorm(abs(Z_join), lower.tail = FALSE)

  result$p_join_vec <- p_join_vec
  result$p_join <- SMiXcan::safe_ACAT(p_join_vec)
  result$mode <- "joint"
  result
}

# Read and standardize the three-cell prediction weights.
weights <- fread(weight_path)
required_weight_cols <- c(
  "gene", "varID", "ref_allele", "eff_allele",
  "weight_cell_1", "weight_cell_2", "weight_cell_3", "type"
)
missing_weight_cols <- setdiff(required_weight_cols, names(weights))
if (length(missing_weight_cols) > 0L) {
  stop(
    "Three-cell weight file is missing columns: ",
    paste(missing_weight_cols, collapse = ", "),
    call. = FALSE
  )
}

weight_columns <- c("weight_cell_1", "weight_cell_2", "weight_cell_3")
nonzero_weight <- rowSums(abs(as.matrix(weights[, ..weight_columns])), na.rm = TRUE) > 0
weights <- weights[nonzero_weight]
weights[, CHR := sub("_.*$", "", varID)]
weights[, POS := suppressWarnings(as.integer(
  tstrsplit(varID, "_", fixed = TRUE, keep = 2L)[[1L]]
))]
weights <- weights[CHR == "chr21" & !is.na(POS)]
if (nrow(weights) == 0L) {
  stop("No nonzero chromosome 21 weights were found.", call. = FALSE)
}

# Read only DRIVE samples that have a non-missing binary phenotype.
phenotype <- fread(phenotype_path)
if (!all(c("subject_index", "affection_status") %in% names(phenotype))) {
  stop("DRIVE phenotype file lacks subject_index or affection_status.", call. = FALSE)
}
phenotype <- phenotype[
  !is.na(subject_index) & affection_status %in% c(0, 1),
  .(subject_index = as.character(subject_index), affection_status)
]
phenotype <- unique(phenotype, by = "subject_index")

genotype_header <- names(fread(genotype_path, nrows = 0L))
annotation_cols <- c("#CHROM", "POS", "ID", "REF", "ALT")
sample_cols <- intersect(phenotype$subject_index, genotype_header)
if (length(sample_cols) == 0L) {
  stop("No DRIVE phenotype IDs match genotype sample columns.", call. = FALSE)
}

genotype <- fread(
  genotype_path,
  select = c(annotation_cols, sample_cols),
  na.strings = c("NA", ".")
)
setnames(genotype, "#CHROM", "CHR")
genotype[, variant_index := .I]

# Match variants and orient every weight to the ALT dosage used in DRIVE.
matched <- merge(
  weights,
  genotype[, .(variant_index, CHR, POS, REF, ALT)],
  by = c("CHR", "POS"),
  allow.cartesian = TRUE
)
matched <- matched[
  (eff_allele == ALT & ref_allele == REF) |
    (eff_allele == REF & ref_allele == ALT)
]
matched[, allele_sign := fifelse(eff_allele == ALT, 1, -1)]
for (column in c("weight_cell_1", "weight_cell_2", "weight_cell_3")) {
  set(matched, j = column, value = matched[[column]] * matched$allele_sign)
}
matched <- unique(matched, by = c("gene", "variant_index"))
if (nrow(matched) == 0L) {
  stop("No allele-aligned DRIVE variants overlap the three-cell weights.", call. = FALSE)
}

# Form the sample-by-variant dosage matrix and mean-impute sporadic missingness.
variant_indices <- sort(unique(matched$variant_index))
dosage <- t(as.matrix(genotype[variant_indices, ..sample_cols]))
storage.mode(dosage) <- "double"
rownames(dosage) <- sample_cols
colnames(dosage) <- as.character(variant_indices)
for (j in seq_len(ncol(dosage))) {
  missing <- is.na(dosage[, j])
  if (any(missing)) {
    dosage[missing, j] <- mean(dosage[, j], na.rm = TRUE)
  }
}
finite_variants <- colSums(is.finite(dosage)) == nrow(dosage)
dosage <- dosage[, finite_variants, drop = FALSE]
variant_indices <- as.integer(colnames(dosage))
matched <- matched[variant_index %in% variant_indices]

outcome <- phenotype[match(rownames(dosage), subject_index), affection_status]
sample_keep <- !is.na(outcome)
dosage <- dosage[sample_keep, , drop = FALSE]
outcome <- outcome[sample_keep]
n0 <- sum(outcome == 0)
n1 <- sum(outcome == 1)
if (n0 == 0L || n1 == 0L) {
  stop("DRIVE analysis requires both cases and controls.", call. = FALSE)
}

message(
  "Running DRIVE SNP GWAS for ", ncol(dosage), " variants in ",
  nrow(dosage), " samples (", n1, " cases; ", n0, " controls)."
)
gwas <- as.data.table(SMiXcan::run_gwas(
  dosage,
  outcome,
  stats::binomial(),
  method_binomial = "glm"
))
gwas[, variant_index := variant_indices]

genes <- unique(matched[, .(gene, type)])
result_rows <- vector("list", nrow(genes))

for (g in seq_len(nrow(genes))) {
  gene_id <- genes$gene[g]
  gene_type <- genes$type[g]
  gene_weights <- matched[gene == gene_id]
  gene_weights <- gene_weights[match(
    intersect(variant_indices, gene_weights$variant_index),
    variant_index
  )]
  gene_gwas <- gwas[match(gene_weights$variant_index, variant_index)]
  valid <- is.finite(gene_gwas$Beta) & is.finite(gene_gwas$se_Beta) &
    gene_gwas$se_Beta > 0
  gene_weights <- gene_weights[valid]
  gene_gwas <- gene_gwas[valid]

  if (nrow(gene_weights) == 0L) {
    next
  }

  x_gene <- dosage[, match(gene_weights$variant_index, variant_indices), drop = FALSE]
  W <- as.matrix(gene_weights[, .(
    weight_cell_1, weight_cell_2, weight_cell_3
  )])
  predicted <- x_gene %*% W

  mixcan <- mixcan_assoc_test_k(outcome, predicted, family = "binomial")
  smixcan <- tryCatch(
    SMiXcan::SMiXcan_assoc_test_K(
      W = W,
      gwas_results = list(Beta = gene_gwas$Beta, se_Beta = gene_gwas$se_Beta),
      x_g = x_gene,
      n0 = n0,
      n1 = n1,
      family = "binomial"
    ),
    error = function(e) NULL
  )
  if (is.null(smixcan)) {
    next
  }

  result_rows[[g]] <- data.table(
    gene_id = gene_id,
    type = gene_type,
    input_snp_num = nrow(gene_weights),
    p_m_sep_1 = mixcan$p_sep[1],
    p_m_sep_2 = mixcan$p_sep[2],
    p_m_sep_3 = mixcan$p_sep[3],
    p_m_sep = mixcan$p_sep_combined,
    p_m_join_1 = mixcan$p_join_vec[1],
    p_m_join_2 = mixcan$p_join_vec[2],
    p_m_join_3 = mixcan$p_join_vec[3],
    p_m_join = mixcan$p_join,
    p_s_sep_1 = smixcan$p_sep[1],
    p_s_sep_2 = smixcan$p_sep[2],
    p_s_sep_3 = smixcan$p_sep[3],
    p_s_sep = SMiXcan::safe_ACAT(smixcan$p_sep),
    p_s_join_1 = smixcan$p_join_vec[1],
    p_s_join_2 = smixcan$p_join_vec[2],
    p_s_join_3 = smixcan$p_join_vec[3],
    p_s_join = smixcan$p_join,
    mixcan_mode = mixcan$mode,
    smixcan_mode = smixcan$mode
  )

  if (g %% 10L == 0L || g == nrow(genes)) {
    message("Processed ", g, " of ", nrow(genes), " DRIVE-overlapping genes.")
  }
}

drive_result <- rbindlist(result_rows, fill = TRUE)
if (nrow(drive_result) == 0L) {
  stop("No three-cell DRIVE gene results were produced.", call. = FALSE)
}
fwrite(drive_result, output_path)

message("Three-cell DRIVE result: ", output_path)
message("Genes written: ", nrow(drive_result))

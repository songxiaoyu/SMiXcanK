# Run S-MiXcan association for HERMES_HF HF using PWAS heart-protein weights.
#
# Inputs from previous steps:
#   2_HERMES_HF_prepare_data_pwas.R:
#     Results/hermes_hf_pwas/hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_hf_input/*.rds
#   3_1000Genome_keep_eur_plink_hermes_hf_pwas.sh:
#     Results/hermes_hf_pwas/hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin/hermes_hf_filtered_id/*.bim/*.raw
#
# HERMES_HF HF is treated as a binary trait by default.
# Provided counts:
#   cases    = 132176
#   controls = 1553537
# This active version uses the package-level association function directly.

library(data.table)
library(SMiXcan)

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hf_pwas",
  "hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
hermes_hf_input_dir <- file.path(workspace_dir, "hermes_hf_input")
ld_input_dir <- file.path(workspace_dir, "hermes_hf_filtered_id")
result_dir <- file.path(workspace_dir, "hermes_hf_result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

family <- "binomial"
n1 <- 132176
n0 <- 1553537
args <- commandArgs(trailingOnly = TRUE)
chr_list <- if (length(args)) as.integer(args) else 1:22
assoc_reg_scale <- as.numeric(Sys.getenv("HERMES_HF_ASSOC_REG_SCALE", unset = "0.2"))
scale_tag <- gsub("\\.", "p", sprintf("%.3f", assoc_reg_scale))
result_tag <- paste0("fixed_", scale_tag)
gene_log_every <- 50
# When a chromosome argument is supplied, this script is being called by the
# parallel wrapper. In that mode each worker writes one chr result; the wrapper
# combines all chromosomes after every worker finishes.
write_combined <- length(args) == 0

cat("Using fixed regularization scale =", assoc_reg_scale, "\n")
cat("Result tag =", result_tag, "\n")

run_gene <- function(selected, X_ref, chr) {
  selected_snp_id <- selected$varID
  available <- selected_snp_id[selected_snp_id %in% colnames(X_ref)]
  selected <- selected[selected$varID %in% available, , drop = FALSE]
  if (!nrow(selected)) {
    return(NULL)
  }

  X_ref_filtered <- X_ref[, selected$varID, drop = FALSE]
  W <- cbind(
    selected$weight_cardiomyocytes,
    selected$weight_other
  )
  gwas_results <- list(
    Beta = selected$beta.Gwas,
    se_Beta = selected$SE.Gwas
  )

  # SMiXcan_assoc_test_K() expects SNP x cell-type weights, GWAS beta/se for
  # the same SNP order, and a 1000 Genomes reference dosage matrix.
  fit <- SMiXcan_assoc_test_K(
    W = W,
    gwas_results = gwas_results,
    x_g = X_ref_filtered,
    n0 = n0,
    n1 = n1,
    family = family,
    regularization = "fixed",
    reg_scale = assoc_reg_scale
  )

  # p_join_pre_reg is a diagnostic ACAT combination of the marginal
  # single-cell-type p-values before the joint regularized inversion.
  pre_reg_p_join <- SMiXcan:::safe_ACAT(fit$p_sep)
  scale_col <- paste0("p_join_reg_scale_", scale_tag)

  out <- data.frame(
    gene_id = selected$gene_id[1],
    gene_name = selected$gene_name[1],
    chr = chr,
    type = selected$type[1],
    input_snp_num = nrow(selected),
    assoc_reg_mode = "fixed",
    assoc_reg_scale_selected = fit$reg_scale_selected,
    assoc_reg_condition = fit$reg_condition,
    Z_cardiomyocytes = fit$Z_join[1],
    p_cardiomyocytes = fit$p_join_vec[1],
    Z_other = fit$Z_join[2],
    p_other = fit$p_join_vec[2],
    p_join = fit$p_join,
    Z_cardiomyocytes_pre_reg = fit$Z_sep[1],
    p_cardiomyocytes_pre_reg = fit$p_sep[1],
    Z_other_pre_reg = fit$Z_sep[2],
    p_other_pre_reg = fit$p_sep[2],
    p_join_pre_reg = pre_reg_p_join,
    Z_cardiomyocytes_reg = fit$Z_join[1],
    p_cardiomyocytes_reg = fit$p_join_vec[1],
    Z_other_reg = fit$Z_join[2],
    p_other_reg = fit$p_join_vec[2],
    p_join_reg = fit$p_join,
    check.names = FALSE
  )
  out[[scale_col]] <- fit$p_join
  out
}

for (chr in chr_list) {
  cat("Processing chromosome", chr, "\n")
  result_path <- file.path(
    result_dir,
    sprintf("hermes_hf_chr%d_result_pwas_%s.csv", chr, result_tag)
  )
  if (file.exists(result_path)) {
    cat("  --> Existing result found, skipping:", result_path, "\n")
    next
  }

  mw_gwas_input_path <- file.path(
    hermes_hf_input_dir,
    sprintf("chr%d_mw_gwas_input_hermes_hf_pwas.rds", chr)
  )
  ld_snp_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_hermes_hf_pwas.bim", chr))
  x_ref_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_012_hermes_hf_pwas.raw", chr))

  if (!file.exists(mw_gwas_input_path) || !file.exists(ld_snp_path) || !file.exists(x_ref_path)) {
    cat("  --> Missing input files, skipping\n")
    next
  }

  mw_gwas_input <- readRDS(mw_gwas_input_path)
  if (!"gene_id" %in% names(mw_gwas_input) && "protein_id" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_id <- mw_gwas_input$protein_id
  }
  if (!"gene_name" %in% names(mw_gwas_input) && "protein_name" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_name <- mw_gwas_input$protein_name
  }
  if (!all(c("gene_id", "gene_name") %in% names(mw_gwas_input))) {
    stop("Step 4 input must contain gene_id/gene_name or protein_id/protein_name columns.", call. = FALSE)
  }
  if (nrow(mw_gwas_input) == 0) {
    cat("  --> Empty merged weight/GWAS input, skipping\n")
    next
  }

  ld_snp <- fread(ld_snp_path, header = FALSE)
  ref_snp_id <- ld_snp$V2

  x_ref_dt <- fread(x_ref_path)
  if (ncol(x_ref_dt) <= 6) {
    cat("  --> Empty LD reference raw file, skipping\n")
    next
  }
  # PLINK2 --export A writes six sample metadata columns followed by SNP
  # dosage columns named like varID_allele. Remove the allele suffix to match
  # the varID used in the weights/GWAS input.
  X_ref <- as.matrix(x_ref_dt[, 7:ncol(x_ref_dt)])
  colnames(X_ref) <- sub("_[^_]+$", "", colnames(X_ref))

  filtered_mw_gwas_input <- mw_gwas_input[mw_gwas_input$varID %in% intersect(ref_snp_id, colnames(X_ref)), ]
  if (nrow(filtered_mw_gwas_input) == 0) {
    cat("  --> No overlap with LD reference SNPs, skipping\n")
    next
  }

  split_df <- split(as.data.frame(filtered_mw_gwas_input), filtered_mw_gwas_input$gene_id)
  result_list <- vector("list", length(split_df))
  for (g in seq_along(split_df)) {
    gene_id <- names(split_df)[g]
    if (gene_log_every > 0 && (g == 1 || g %% gene_log_every == 0 || g == length(split_df))) {
      cat("  Gene", g, "of", length(split_df), ":", gene_id, "\n")
    }
    result_list[[g]] <- run_gene(split_df[[g]], X_ref, chr)
  }

  real_result <- rbindlist(result_list, fill = TRUE)
  write.csv(real_result, result_path, row.names = FALSE)
}

available_results <- file.path(
  result_dir,
  sprintf("hermes_hf_chr%d_result_pwas_%s.csv", chr_list, result_tag)
)
available_results <- available_results[file.exists(available_results)]
if (write_combined && length(available_results)) {
  combined <- rbindlist(lapply(available_results, fread), fill = TRUE)
  write.csv(
    combined,
    file.path(result_dir, sprintf("hermes_hf_result_pwas_%s.csv", result_tag)),
    row.names = FALSE
  )
}

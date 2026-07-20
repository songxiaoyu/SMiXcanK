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

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
analysis_dir <- file.path(paper_dir, "Results", "2pi_workspace")
bcac_input_dir <- file.path(analysis_dir, "bcac2020_input")
ld_input_dir <- file.path(analysis_dir, "bcac2020_filtered_id")
result_dir <- file.path(analysis_dir, "bcac2020_result")

# STEP1: load data--------
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

n1 <- 133384

n0 <- 113789 + 18908


# Set chromosome
for(chr in 1:22){
  print('CHR')
  print(chr)
  # Load pre-merged data
  mw_gwas_input_path <- file.path(bcac_input_dir, sprintf("chr%d_mw_gwas_input_bcac2020_pi2.rds", chr))
  mw_gwas_input <- readRDS(mw_gwas_input_path)

  # Load LD matrix, ref genome, SNP info
  LD_snp_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_pi2.bim", chr))
  X_ref_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_012_pi2.raw", chr))

  ld_snp <- fread(LD_snp_path, header = FALSE)
  ref_snp_id <- ld_snp$V2

  # Read genotype matrix and set proper colnames
  X_ref_dt <- fread(X_ref_path)
  X_ref <- as.matrix(X_ref_dt[, 7:ncol(X_ref_dt)])
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
    W <- cbind(W1, W2)
    filtered_list[[gene]] <- list(W = W, selected_snp = gene_df)
  }


  # STEP2 Run S-MiXcan----
  G <- length(filtered_list)
  real_result = data.frame(matrix(ncol = 9, nrow = G))
  colnames(real_result) <- c('gene_id','chr','type','input_snp_num','Z_1','p_1','Z_2','p_2','p_join')
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
    real_result[g, c('gene_id','chr', 'type', 'input_snp_num')] =c(filtered_list[[gene]]$selected_snp[1,c('gene','CHR','type')], nrow(W))
    real_result[g, c('Z_1','p_1','Z_2','p_2','p_join')] <- c(c(rbind(S_MiXcan_results$Z_join, S_MiXcan_results$p_join_vec)), S_MiXcan_results$p_join)
  }

  result_path <- file.path(result_dir, sprintf("bcac2020_chr%d_result_pi2.csv", chr))
  write.csv(real_result, result_path, row.names = FALSE)
}

# Combine 22 chr and save the result
combined <- do.call(rbind, lapply(1:22, function(chr) {
  result_path <- file.path(result_dir, sprintf("bcac2020_chr%d_result_pi2.csv", chr))
  read.csv(result_path)
}))
combined_path <- file.path(result_dir, "bcac2020_result_pi2.csv")
write.csv(combined, combined_path, row.names = FALSE)

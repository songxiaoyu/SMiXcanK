# Test fixed regularization scales for HERMES HF PWAS.
#
# This diagnostic reuses existing step 2/3 HF inputs and does not change the
# main HF workflow defaults.
#
# Optional output tag:
#   HERMES_HF_REG_SENS_TAG=my_tag Rscript Analysis_code/HERMES_HF/7_optional_test_regularization_scales_pwas.R
#   Rscript Analysis_code/HERMES_HF/7_optional_test_regularization_scales_pwas.R my_tag
#
# With no tag, the historical output file names are unchanged.

library(data.table)
library(bacon)
library(SMiXcan)

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hf_pwas",
  "hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
input_dir <- file.path(workspace_dir, "hermes_hf_input")
ld_input_dir <- file.path(workspace_dir, "hermes_hf_filtered_id")
out_dir <- file.path(workspace_dir, "hermes_hf_result", "regularization_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
output_tag <- Sys.getenv("HERMES_HF_REG_SENS_TAG", unset = "")
if (!nzchar(output_tag) && length(args) >= 1) {
  output_tag <- args[1]
}
if (nzchar(output_tag)) {
  output_tag <- gsub("[^A-Za-z0-9_.-]+", "_", output_tag)
  output_suffix <- paste0("_", output_tag)
} else {
  output_suffix <- ""
}

n1 <- 132176
n0 <- 1553537
family <- "binomial"
chr_list <- 1:22
scales <- c(0.001, 0.002, 0.005, 0.01, 0.03, 0.05, 0.1)

run_gene_scale <- function(selected, X_ref, chr, scale) {
  selected <- selected[selected$varID %in% colnames(X_ref), , drop = FALSE]
  if (!nrow(selected)) {
    return(NULL)
  }

  X_ref_filtered <- X_ref[, selected$varID, drop = FALSE]
  W <- cbind(selected$weight_cardiomyocytes, selected$weight_other)
  gwas_results <- list(Beta = selected$beta.Gwas, se_Beta = selected$SE.Gwas)

  fit <- SMiXcan_assoc_test_K(
    W = W,
    gwas_results = gwas_results,
    x_g = X_ref_filtered,
    n0 = n0,
    n1 = n1,
    family = family,
    regularization = "fixed",
    reg_scale = scale
  )

  data.table(
    gene_id = selected$gene_id[1],
    gene_name = selected$gene_name[1],
    chr = chr,
    type = selected$type[1],
    input_snp_num = nrow(selected),
    reg_scale = scale,
    p_join = fit$p_join,
    p_cardiomyocytes = fit$p_join_vec[1],
    p_other = fit$p_join_vec[2],
    p_join_pre_reg = SMiXcan:::safe_ACAT(fit$p_sep)
  )
}

run_chr <- function(chr) {
  cat("Processing chr", chr, "\n")
  mw_gwas_input_path <- file.path(
    input_dir,
    sprintf("chr%d_mw_gwas_input_hermes_hf_pwas.rds", chr)
  )
  ld_snp_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_hermes_hf_pwas.bim", chr))
  x_ref_path <- file.path(ld_input_dir, sprintf("filtered_chr%d_hg38_012_hermes_hf_pwas.raw", chr))

  if (!file.exists(mw_gwas_input_path) || !file.exists(ld_snp_path) || !file.exists(x_ref_path)) {
    warning("Missing input for chr", chr)
    return(NULL)
  }

  mw_gwas_input <- readRDS(mw_gwas_input_path)
  if (!"gene_id" %in% names(mw_gwas_input) && "protein_id" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_id <- mw_gwas_input$protein_id
  }
  if (!"gene_name" %in% names(mw_gwas_input) && "protein_name" %in% names(mw_gwas_input)) {
    mw_gwas_input$gene_name <- mw_gwas_input$protein_name
  }

  ld_snp <- fread(ld_snp_path, header = FALSE)
  ref_snp_id <- ld_snp$V2

  x_ref_dt <- fread(x_ref_path)
  X_ref <- as.matrix(x_ref_dt[, 7:ncol(x_ref_dt)])
  colnames(X_ref) <- sub("_[^_]+$", "", colnames(X_ref))

  filtered <- mw_gwas_input[mw_gwas_input$varID %in% intersect(ref_snp_id, colnames(X_ref)), ]
  if (!nrow(filtered)) {
    return(NULL)
  }

  split_df <- split(as.data.frame(filtered), filtered$gene_id)
  chr_result <- rbindlist(lapply(split_df, function(selected) {
    rbindlist(lapply(scales, function(scale) {
      run_gene_scale(selected, X_ref, chr, scale)
    }), fill = TRUE)
  }), fill = TRUE)

  chr_path <- file.path(
    out_dir,
    sprintf("hermes_hf_chr%d_regularization_sensitivity%s.csv", chr, output_suffix)
  )
  fwrite(chr_result, chr_path)
  chr_result
}

all_result <- rbindlist(lapply(chr_list, run_chr), fill = TRUE)
all_path <- file.path(
  out_dir,
  sprintf("hermes_hf_regularization_sensitivity_all%s.csv", output_suffix)
)
fwrite(all_result, all_path)

calc_lambda_gc <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p < 1]
  if (!length(p)) {
    return(NA_real_)
  }
  median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

calc_bacon_inflation <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p < 1]
  if (!length(p)) {
    return(NA_real_)
  }
  z <- qnorm(1 - pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps),
             lower.tail = FALSE)
  as.numeric(inflation(bacon(z, na.exclude = TRUE)))
}

calc_bacon_bias <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p < 1]
  if (!length(p)) {
    return(NA_real_)
  }
  z <- qnorm(1 - pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps),
             lower.tail = FALSE)
  as.numeric(bias(bacon(z, na.exclude = TRUE)))
}

summary_result <- all_result[
  ,
  .(
    n_tested = .N,
    min_p = min(p_join, na.rm = TRUE),
    lambda = calc_bacon_inflation(p_join),
    lambda_gc = calc_lambda_gc(p_join),
    bacon_bias = calc_bacon_bias(p_join),
    fdr_lt_0p05 = sum(p.adjust(p_join, method = "fdr") < 0.05, na.rm = TRUE),
    fdr_lt_0p10 = sum(p.adjust(p_join, method = "fdr") < 0.10, na.rm = TRUE),
    nominal_p_lt_0p05 = sum(p_join < 0.05, na.rm = TRUE),
    median_p = median(p_join, na.rm = TRUE)
  ),
  by = reg_scale
][order(reg_scale)]

summary_path <- file.path(
  out_dir,
  sprintf("hermes_hf_regularization_sensitivity_summary%s.csv", output_suffix)
)
fwrite(summary_result, summary_path)

cat("Output tag:", ifelse(nzchar(output_tag), output_tag, "<none>"), "\n")
cat("Wrote full result:", all_path, "\n")
cat("Wrote summary:", summary_path, "\n")
print(summary_result)

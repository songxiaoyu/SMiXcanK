#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

project_dir <- get_env("PAPER_SMIXCAN_DIR",
                       "/Users/admin/Library/CloudStorage/Dropbox/Paper_SMiXcan")
output_prefix <- get_env(
  "PWAS_OUTPUT_PREFIX",
  "_moderate_100kb_r2_0.99_alpha0.5_lambdamin_predixcan_with_intercepts"
)
dosage_dir <- get_env(
  "PWAS_DOSAGE_DIR",
  file.path(project_dir,
            "New generated files/codes/by_chr_heart_left_ventricle_351_samples_variant_id")
)
weights_file <- get_env(
  "PWAS_PREDIXCAN_WEIGHTS_FILE",
  file.path(project_dir,
            "Results/heart_protein_weights/training_model_weights",
            paste0("predixcan_tissue_weights_heart_protein", output_prefix, ".csv"))
)
intercepts_file <- get_env(
  "PWAS_INTERCEPTS_FILE",
  file.path(project_dir,
            "Results/heart_protein_weights/training_model_weights",
            paste0("intercepts_heart_protein_cardiomyocytes_other",
                   output_prefix, ".csv"))
)
out_dir <- get_env(
  "PWAS_PREDIXCAN_OUT_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/predicted_protein_HLV351_predixcan_with_intercepts")
)
chromosomes <- as.integer(strsplit(get_env("PWAS_CHR", paste(1:22, collapse = ",")),
                                   ",", fixed = TRUE)[[1]])
chromosomes <- chromosomes[!is.na(chromosomes)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Dosage dir: ", dosage_dir)
message("PrediXcan tissue weights: ", weights_file)
message("Intercepts: ", intercepts_file)
message("Output dir: ", out_dir)

weights <- fread(weights_file)
intercepts <- fread(intercepts_file)

required_weight_cols <- c("gene_id", "gene_name", "varID", "chr",
                          "dosed_allele", "weight_tissue")
missing_weight_cols <- setdiff(required_weight_cols, names(weights))
if (length(missing_weight_cols)) {
  stop("PrediXcan weights file is missing columns: ",
       paste(missing_weight_cols, collapse = ", "), call. = FALSE)
}

required_intercept_cols <- c("gene_id", "gene_name", "chr", "intercept_tissue")
missing_intercept_cols <- setdiff(required_intercept_cols, names(intercepts))
if (length(missing_intercept_cols)) {
  stop("Intercepts file is missing columns: ",
       paste(missing_intercept_cols, collapse = ", "), call. = FALSE)
}

prediction_files <- character()
match_reports <- list()

for (chr in chromosomes) {
  chr_value <- chr
  message("===== Predicting PrediXcan chr", chr, " =====")

  dosage_file <- file.path(
    dosage_dir,
    paste0("chr", chr, "_HLV351_dosage_nomiss_variant_id.raw")
  )
  if (!file.exists(dosage_file)) {
    dosage_file <- file.path(dosage_dir, paste0("chr", chr, "_HLV351_dosage_nomiss.raw"))
  }
  if (!file.exists(dosage_file)) {
    stop("Missing dosage file for chr", chr, ": ", dosage_file, call. = FALSE)
  }

  intercept_chr <- intercepts[as.integer(chr) == chr_value]
  weights_chr <- weights[as.integer(chr) == chr_value]
  if (!nrow(intercept_chr)) {
    warning("No intercept rows for chr", chr, "; skipping")
    next
  }

  raw_header <- names(fread(dosage_file, nrows = 0))
  id_cols <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), raw_header)

  weights_chr[, dosage_col := paste0(varID, "_", dosed_allele)]
  weights_chr[!(dosage_col %in% raw_header), dosage_col := NA_character_]

  match_report_chr <- weights_chr[, .(
    chr = chr_value,
    n_weight_rows = .N,
    n_with_dosage_col = sum(!is.na(dosage_col)),
    n_missing_dosage_col = sum(is.na(dosage_col))
  )]
  match_reports[[as.character(chr)]] <- match_report_chr

  missing_rows <- weights_chr[is.na(dosage_col),
                              .(gene_id, gene_name, varID, chr,
                                dosed_allele, weight_tissue)]
  if (nrow(missing_rows)) {
    fwrite(missing_rows,
           file.path(out_dir, paste0("missing_predixcan_weight_snps_chr",
                                     chr, ".csv")))
  }

  use_weights <- weights_chr[!is.na(dosage_col)]
  dosage_cols <- unique(use_weights$dosage_col)
  dosage <- fread(dosage_file, select = unique(c(id_cols, dosage_cols)))
  sample_dt <- dosage[, .(FID, IID)]
  sample_dt[, sample_id := IID]
  n_samples <- nrow(sample_dt)

  if (length(dosage_cols)) {
    dosage_mat <- as.matrix(dosage[, ..dosage_cols])
    storage.mode(dosage_mat) <- "numeric"
    colnames(dosage_mat) <- dosage_cols
  } else {
    dosage_mat <- matrix(numeric(0), nrow = n_samples, ncol = 0)
  }

  pred_list <- vector("list", nrow(intercept_chr))
  for (i in seq_len(nrow(intercept_chr))) {
    gene <- intercept_chr$gene_id[i]
    gene_weights <- use_weights[gene_id == gene]
    pred_tissue <- rep(intercept_chr$intercept_tissue[i], n_samples)

    if (nrow(gene_weights)) {
      x <- dosage_mat[, gene_weights$dosage_col, drop = FALSE]
      pred_tissue <- pred_tissue + as.numeric(x %*% gene_weights$weight_tissue)
    }

    pred_list[[i]] <- data.table(
      FID = sample_dt$FID,
      IID = sample_dt$IID,
      sample_id = sample_dt$sample_id,
      gene_id = gene,
      gene_name = intercept_chr$gene_name[i],
      chr = chr_value,
      pred_tissue_predixcan = pred_tissue
    )
  }

  pred_chr <- rbindlist(pred_list, use.names = TRUE)
  pred_file <- file.path(out_dir, paste0("predicted_heart_protein_HLV351_predixcan_chr",
                                         chr, "_with_intercepts.csv"))
  fwrite(pred_chr, pred_file)
  prediction_files <- c(prediction_files, pred_file)

  message("Rows written chr", chr, ": ", nrow(pred_chr))
  message("Matched PrediXcan weight rows chr", chr, ": ",
          match_report_chr$n_with_dosage_col, " / ", match_report_chr$n_weight_rows)

  rm(dosage, dosage_mat, pred_chr, pred_list)
  gc()
}

match_report <- rbindlist(match_reports, use.names = TRUE, fill = TRUE)
fwrite(match_report, file.path(out_dir, "predixcan_snp_matching_report_by_chr.csv"))

message("Combining PrediXcan predictions...")
pred_all <- rbindlist(lapply(prediction_files, fread), use.names = TRUE)
combined_file <- file.path(out_dir,
                           "predicted_heart_protein_HLV351_predixcan_all_chr_with_intercepts.csv")
fwrite(pred_all, combined_file)

summary_file <- file.path(out_dir, "predixcan_prediction_summary.txt")
summary_lines <- c(
  paste0("PrediXcan prediction date: ", Sys.time()),
  paste0("Dosage dir: ", dosage_dir),
  paste0("PrediXcan weights file: ", weights_file),
  paste0("Intercepts file: ", intercepts_file),
  paste0("Combined prediction file: ", combined_file),
  paste0("Samples: ", uniqueN(pred_all$sample_id)),
  paste0("Genes: ", uniqueN(pred_all$gene_id)),
  paste0("Rows: ", nrow(pred_all)),
  paste0("Total PrediXcan weight rows: ", nrow(weights)),
  paste0("Matched PrediXcan weight rows: ", sum(match_report$n_with_dosage_col)),
  paste0("Missing dosage-column PrediXcan weight rows: ",
         sum(match_report$n_missing_dosage_col)),
  "Formula: pred_tissue_predixcan = intercept_tissue + sum(dosage * weight_tissue)"
)
writeLines(summary_lines, summary_file)

message("Saved combined PrediXcan predictions: ", combined_file)
message("Saved matching report: ",
        file.path(out_dir, "predixcan_snp_matching_report_by_chr.csv"))
message("Saved summary: ", summary_file)

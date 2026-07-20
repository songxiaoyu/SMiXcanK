#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

donor_id <- function(x) {
  vapply(strsplit(x, "-", fixed = TRUE), function(parts) {
    paste(parts[seq_len(min(2L, length(parts)))], collapse = "-")
  }, character(1))
}

project_dir <- get_env("PAPER_SMIXCAN_DIR",
                       "/Users/admin/Library/CloudStorage/Dropbox/Paper_SMiXcan")
prediction_dir <- get_env(
  "PWAS_PREDICTION_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts")
)
cell_prediction_file <- get_env(
  "PWAS_CELL_PREDICTION_FILE",
  file.path(prediction_dir,
            "predicted_heart_protein_HLV351_all_chr_with_intercepts.csv")
)
pi_file <- get_env(
  "PWAS_PI_FILE",
  file.path(project_dir, "Heart/GTEx_Pi_Estimate/RNA_all_2celltypes_351_pi.tsv")
)
out_file <- get_env(
  "PWAS_TISSUE_PREDICTION_FILE",
  file.path(prediction_dir,
            "predicted_heart_protein_HLV351_tissue_level_with_intercepts.csv")
)
summary_file <- get_env(
  "PWAS_TISSUE_PREDICTION_SUMMARY",
  file.path(prediction_dir, "tissue_level_prediction_summary.txt")
)

message("Cell-level prediction file: ", cell_prediction_file)
message("PI file: ", pi_file)
message("Output file: ", out_file)

pred <- fread(cell_prediction_file)
pi <- fread(pi_file)

if (!all(c("pred_cardiomyocytes", "pred_other") %in% names(pred))) {
  stop("Prediction file must contain pred_cardiomyocytes and pred_other.",
       call. = FALSE)
}
if (!all(c("Cardiomyocytes", "Other") %in% names(pi))) {
  stop("PI file must contain Cardiomyocytes and Other columns.", call. = FALSE)
}

pi_id_col <- setdiff(names(pi), c("Cardiomyocytes", "Other"))[1]
if (is.na(pi_id_col)) {
  stop("Could not identify donor ID column in PI file.", call. = FALSE)
}
setnames(pi, pi_id_col, "donor_id")
pi[, donor_id := as.character(donor_id)]

if (!"sample_id" %in% names(pred)) {
  if ("IID" %in% names(pred)) {
    pred[, sample_id := IID]
  } else {
    stop("Prediction file must contain sample_id or IID.", call. = FALSE)
  }
}
pred[, donor_id := donor_id(sample_id)]

missing_pi <- setdiff(unique(pred$donor_id), unique(pi$donor_id))
if (length(missing_pi)) {
  stop("Missing PI rows for ", length(missing_pi), " donors. Examples: ",
       paste(head(missing_pi, 10), collapse = ", "), call. = FALSE)
}

pi[, pi_sum := Cardiomyocytes + Other]
bad_pi <- pi[!is.finite(pi_sum) | abs(pi_sum - 1) > 1e-4]
if (nrow(bad_pi)) {
  warning("Some PI rows do not sum to 1 within 1e-4. Max deviation: ",
          max(abs(bad_pi$pi_sum - 1), na.rm = TRUE))
}
pi[, pi_sum := NULL]

pred <- merge(
  pred,
  pi[, .(donor_id, pi_cardiomyocytes = Cardiomyocytes, pi_other = Other)],
  by = "donor_id",
  all.x = TRUE,
  sort = FALSE
)

pred[, pred_tissue := pi_cardiomyocytes * pred_cardiomyocytes +
       pi_other * pred_other]

out <- pred[, .(
  FID, IID, sample_id, donor_id, gene_id, gene_name, chr,
  pi_cardiomyocytes, pi_other,
  pred_cardiomyocytes, pred_other, pred_tissue
)]

fwrite(out, out_file)

summary_lines <- c(
  paste0("Tissue-level prediction date: ", Sys.time()),
  paste0("Cell-level prediction file: ", cell_prediction_file),
  paste0("PI file: ", pi_file),
  paste0("Output file: ", out_file),
  paste0("Samples: ", uniqueN(out$sample_id)),
  paste0("Donors: ", uniqueN(out$donor_id)),
  paste0("Genes: ", uniqueN(out$gene_id)),
  paste0("Rows: ", nrow(out)),
  paste0("Formula: pred_tissue = pi_cardiomyocytes * pred_cardiomyocytes + pi_other * pred_other")
)
writeLines(summary_lines, summary_file)

message("Saved tissue-level predictions: ", out_file)
message("Saved summary: ", summary_file)

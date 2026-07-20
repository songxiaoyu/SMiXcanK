#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

donor_id <- function(x) {
  x <- as.character(x)
  vapply(strsplit(x, "-", fixed = TRUE), function(parts) {
    paste(parts[seq_len(min(2L, length(parts)))], collapse = "-")
  }, character(1))
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  x <- x[keep]
  y <- y[keep]
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  stats::cor(x, y, method = method)
}

project_dir <- get_env("PAPER_SMIXCAN_DIR",
                       "/Users/admin/Library/CloudStorage/Dropbox/Paper_SMiXcan")
prediction_dir_smixcan <- get_env(
  "PWAS_SMIXCAN_PRED_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts")
)
prediction_dir_predixcan <- get_env(
  "PWAS_PREDIXCAN_PRED_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/predicted_protein_HLV351_predixcan_with_intercepts")
)
smixcan_file <- get_env(
  "PWAS_SMIXCAN_TISSUE_PRED_FILE",
  file.path(prediction_dir_smixcan,
            "predicted_heart_protein_HLV351_tissue_level_with_intercepts.csv")
)
predixcan_file <- get_env(
  "PWAS_PREDIXCAN_PRED_FILE",
  file.path(prediction_dir_predixcan,
            "predicted_heart_protein_HLV351_predixcan_all_chr_with_intercepts.csv")
)
expression_file <- get_env(
  "PWAS_EXPRESSION_GCT",
  file.path(project_dir, "Heart/Data/gene_tpm_2022-06-06_v10_heart_left_ventricle.gct")
)
out_dir <- get_env(
  "PWAS_COR_OUT_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/prediction_expression_correlation_HLV351")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("SMiXcan prediction: ", smixcan_file)
message("PrediXcan prediction: ", predixcan_file)
message("Expression GCT: ", expression_file)
message("Output dir: ", out_dir)

smixcan <- fread(smixcan_file, select = c("sample_id", "donor_id", "gene_id",
                                          "gene_name", "pred_tissue"))
predixcan <- fread(predixcan_file, select = c("sample_id", "gene_id",
                                              "gene_name", "pred_tissue_predixcan"))
predixcan[, donor_id := donor_id(sample_id)]

smixcan <- smixcan[, .(
  donor_id, gene_id, gene_name,
  method = "SMiXcan",
  predicted = pred_tissue
)]
predixcan <- predixcan[, .(
  donor_id, gene_id, gene_name,
  method = "PrediXcan",
  predicted = pred_tissue_predixcan
)]
pred <- rbindlist(list(smixcan, predixcan), use.names = TRUE)
pred <- pred[!is.na(donor_id) & !is.na(gene_id)]

target_genes <- unique(pred$gene_id)
target_donors <- unique(pred$donor_id)

expr <- fread(expression_file, skip = 2, check.names = FALSE)
expr <- expr[Name %in% target_genes]
sample_cols <- setdiff(names(expr), c("Name", "Description"))
sample_donors <- donor_id(sample_cols)
keep_sample_cols <- sample_cols[sample_donors %in% target_donors]

if (!length(keep_sample_cols)) {
  stop("No expression samples overlap prediction donors.", call. = FALSE)
}

expr_long <- melt(
  expr[, c("Name", "Description", keep_sample_cols), with = FALSE],
  id.vars = c("Name", "Description"),
  variable.name = "expression_sample_id",
  value.name = "tpm"
)
setnames(expr_long, c("Name", "Description"), c("gene_id", "expression_gene_name"))
expr_long[, donor_id := donor_id(expression_sample_id)]
expr_long[, tpm := as.numeric(tpm)]

# If a donor has more than one RNA sample, average expression at donor level.
expr_donor <- expr_long[, .(
  expression_gene_name = expression_gene_name[1],
  tpm = mean(tpm, na.rm = TRUE)
), by = .(gene_id, donor_id)]
expr_donor[, log2_tpm_plus1 := log2(tpm + 1)]

dat <- merge(
  pred,
  expr_donor,
  by = c("gene_id", "donor_id"),
  all = FALSE,
  sort = FALSE
)

cor_by_gene <- dat[, .(
  n = .N,
  pearson_tpm = safe_cor(predicted, tpm, "pearson"),
  spearman_tpm = safe_cor(predicted, tpm, "spearman"),
  pearson_log2_tpm_plus1 = safe_cor(predicted, log2_tpm_plus1, "pearson"),
  spearman_log2_tpm_plus1 = safe_cor(predicted, log2_tpm_plus1, "spearman")
), by = .(method, gene_id, gene_name)]

wide <- dcast(
  cor_by_gene,
  gene_id + gene_name ~ method,
  value.var = c("n", "pearson_tpm", "spearman_tpm",
                "pearson_log2_tpm_plus1", "spearman_log2_tpm_plus1")
)
if (all(c("pearson_log2_tpm_plus1_SMiXcan",
          "pearson_log2_tpm_plus1_PrediXcan") %in% names(wide))) {
  wide[, delta_pearson_log2_tpm_plus1 :=
         pearson_log2_tpm_plus1_SMiXcan -
         pearson_log2_tpm_plus1_PrediXcan]
}

cor_file <- file.path(out_dir, "prediction_expression_correlation_by_gene.csv")
wide_file <- file.path(out_dir, "prediction_expression_correlation_by_gene_wide.csv")
plot_file <- file.path(out_dir, "prediction_expression_correlation_boxplot.png")
summary_file <- file.path(out_dir, "prediction_expression_correlation_summary.txt")

fwrite(cor_by_gene, cor_file)
fwrite(wide, wide_file)

plot_dt <- cor_by_gene[!is.na(pearson_log2_tpm_plus1)]
png(plot_file, width = 1800, height = 1400, res = 180)
boxplot(
  pearson_log2_tpm_plus1 ~ method,
  data = plot_dt,
  col = c("#4C78A8", "#F58518"),
  ylab = "Pearson correlation with log2(TPM + 1)",
  xlab = "",
  main = "Predicted Protein vs Real Gene Expression"
)
stripchart(
  pearson_log2_tpm_plus1 ~ method,
  data = plot_dt,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  col = rgb(0, 0, 0, 0.18),
  add = TRUE
)
abline(h = 0, lty = 2, col = "gray40")
dev.off()

summary_dt <- cor_by_gene[, .(
  genes = .N,
  genes_with_nonmissing_pearson_log = sum(!is.na(pearson_log2_tpm_plus1)),
  median_pearson_log2_tpm_plus1 = median(pearson_log2_tpm_plus1, na.rm = TRUE),
  mean_pearson_log2_tpm_plus1 = mean(pearson_log2_tpm_plus1, na.rm = TRUE),
  median_spearman_log2_tpm_plus1 = median(spearman_log2_tpm_plus1, na.rm = TRUE),
  mean_spearman_log2_tpm_plus1 = mean(spearman_log2_tpm_plus1, na.rm = TRUE),
  median_pearson_tpm = median(pearson_tpm, na.rm = TRUE),
  mean_pearson_tpm = mean(pearson_tpm, na.rm = TRUE)
), by = method]

if ("delta_pearson_log2_tpm_plus1" %in% names(wide)) {
  paired <- wide[!is.na(delta_pearson_log2_tpm_plus1)]
  delta_lines <- c(
    paste0("Paired genes: ", nrow(paired)),
    paste0("Median SMiXcan - PrediXcan Pearson log2(TPM+1): ",
           median(paired$delta_pearson_log2_tpm_plus1, na.rm = TRUE)),
    paste0("Mean SMiXcan - PrediXcan Pearson log2(TPM+1): ",
           mean(paired$delta_pearson_log2_tpm_plus1, na.rm = TRUE)),
    paste0("Genes where SMiXcan > PrediXcan: ",
           sum(paired$delta_pearson_log2_tpm_plus1 > 0, na.rm = TRUE)),
    paste0("Genes where PrediXcan > SMiXcan: ",
           sum(paired$delta_pearson_log2_tpm_plus1 < 0, na.rm = TRUE))
  )
} else {
  delta_lines <- "Paired comparison unavailable."
}

summary_lines <- c(
  paste0("Correlation analysis date: ", Sys.time()),
  paste0("SMiXcan prediction file: ", smixcan_file),
  paste0("PrediXcan prediction file: ", predixcan_file),
  paste0("Expression file: ", expression_file),
  paste0("Overlapping donors: ", uniqueN(dat$donor_id)),
  paste0("Overlapping genes: ", uniqueN(dat$gene_id)),
  paste0("Rows used after merge: ", nrow(dat)),
  "",
  "Method summaries:",
  capture.output(print(summary_dt)),
  "",
  "Paired comparison:",
  delta_lines,
  "",
  paste0("Saved per-gene correlations: ", cor_file),
  paste0("Saved wide comparison table: ", wide_file),
  paste0("Saved boxplot: ", plot_file)
)
writeLines(summary_lines, summary_file)

message("Saved per-gene correlations: ", cor_file)
message("Saved wide comparison table: ", wide_file)
message("Saved boxplot: ", plot_file)
message("Saved summary: ", summary_file)

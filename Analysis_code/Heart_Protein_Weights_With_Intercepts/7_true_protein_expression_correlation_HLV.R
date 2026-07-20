#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

normalize_sample_id <- function(x) {
  gsub("\\.", "-", as.character(x))
}

donor_id <- function(x) {
  x <- normalize_sample_id(x)
  vapply(strsplit(x, "-", fixed = TRUE), function(parts) {
    paste(parts[seq_len(min(2L, length(parts)))], collapse = "-")
  }, character(1))
}

safe_cor <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  x <- x[keep]
  y <- y[keep]
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = method)
}

project_dir <- get_env("PAPER_SMIXCAN_DIR",
                       "/Users/admin/Library/CloudStorage/Dropbox/Paper_SMiXcan")
protein_file <- get_env(
  "PWAS_TRUE_PROTEIN_RDATA",
  file.path(project_dir,
            "Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData")
)
expression_file <- get_env(
  "PWAS_EXPRESSION_GCT",
  file.path(project_dir, "Heart/Data/gene_tpm_2022-06-06_v10_heart_left_ventricle.gct")
)
out_dir <- get_env(
  "PWAS_TRUE_COR_OUT_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/true_protein_expression_correlation_HLV")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("True protein RData: ", protein_file)
message("Expression GCT: ", expression_file)
message("Output dir: ", out_dir)

load(protein_file)
if (!exists("Imputed_protein")) {
  stop("Expected object Imputed_protein not found in protein RData.", call. = FALSE)
}
protein <- as.data.table(Imputed_protein)
ann_cols <- c("gene_name", "gene_id", "chr", "start", "end")
protein_sample_cols <- setdiff(names(protein), ann_cols)

protein_long <- melt(
  protein,
  id.vars = intersect(ann_cols, names(protein)),
  measure.vars = protein_sample_cols,
  variable.name = "protein_sample_id",
  value.name = "protein"
)
protein_long[, protein_sample_id := normalize_sample_id(protein_sample_id)]
protein_long[, donor_id := donor_id(protein_sample_id)]
protein_long[, protein := as.numeric(protein)]
protein_donor <- protein_long[, .(
  gene_name = gene_name[1],
  chr = chr[1],
  protein = mean(protein, na.rm = TRUE),
  protein_sample_ids = paste(unique(protein_sample_id), collapse = ";")
), by = .(gene_id, donor_id)]

target_genes <- unique(protein_donor$gene_id)
target_donors <- unique(protein_donor$donor_id)

expr <- fread(expression_file, skip = 2, check.names = FALSE)
expr <- expr[Name %in% target_genes]
sample_cols <- setdiff(names(expr), c("Name", "Description"))
sample_donors <- donor_id(sample_cols)
keep_sample_cols <- sample_cols[sample_donors %in% target_donors]
if (!length(keep_sample_cols)) {
  stop("No RNA expression samples overlap protein donors.", call. = FALSE)
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
expr_donor <- expr_long[, .(
  expression_gene_name = expression_gene_name[1],
  tpm = mean(tpm, na.rm = TRUE),
  expression_sample_ids = paste(unique(expression_sample_id), collapse = ";")
), by = .(gene_id, donor_id)]
expr_donor[, log2_tpm_plus1 := log2(tpm + 1)]

overlap_samples <- merge(
  unique(protein_donor[, .(donor_id, protein_sample_ids)]),
  unique(expr_donor[, .(donor_id, expression_sample_ids)]),
  by = "donor_id",
  all = FALSE
)

dat <- merge(
  protein_donor,
  expr_donor,
  by = c("gene_id", "donor_id"),
  all = FALSE,
  sort = FALSE
)

cor_by_gene <- dat[, .(
  n = .N,
  pearson_tpm = safe_cor(protein, tpm, "pearson"),
  spearman_tpm = safe_cor(protein, tpm, "spearman"),
  pearson_log2_tpm_plus1 = safe_cor(protein, log2_tpm_plus1, "pearson"),
  spearman_log2_tpm_plus1 = safe_cor(protein, log2_tpm_plus1, "spearman")
), by = .(gene_id, gene_name)]

sample_file <- file.path(out_dir, "true_protein_expression_overlap_samples.csv")
merged_file <- file.path(out_dir, "true_protein_expression_merged_values.csv")
cor_file <- file.path(out_dir, "true_protein_expression_correlation_by_gene.csv")
plot_file <- file.path(out_dir, "true_protein_expression_correlation_boxplot.png")
summary_file <- file.path(out_dir, "true_protein_expression_correlation_summary.txt")

fwrite(overlap_samples, sample_file)
fwrite(dat, merged_file)
fwrite(cor_by_gene, cor_file)

png(plot_file, width = 1600, height = 1400, res = 180)
boxplot(
  cor_by_gene$pearson_log2_tpm_plus1,
  col = "#4C78A8",
  ylab = "Pearson correlation with log2(TPM + 1)",
  main = "True Protein vs True RNA Expression"
)
stripchart(
  cor_by_gene$pearson_log2_tpm_plus1,
  vertical = TRUE,
  method = "jitter",
  pch = 16,
  col = rgb(0, 0, 0, 0.2),
  add = TRUE
)
abline(h = 0, lty = 2, col = "gray40")
dev.off()

summary_lines <- c(
  paste0("True protein-expression correlation date: ", Sys.time()),
  paste0("Protein file: ", protein_file),
  paste0("Expression file: ", expression_file),
  paste0("Protein samples: ", length(protein_sample_cols)),
  paste0("Expression samples in GCT: ", length(sample_cols)),
  paste0("Overlapping donors with both protein and RNA: ", nrow(overlap_samples)),
  paste0("Overlapping genes with both protein and RNA: ", uniqueN(dat$gene_id)),
  paste0("Merged gene-donor rows: ", nrow(dat)),
  paste0("Genes with nonmissing Pearson log2(TPM+1): ",
         sum(!is.na(cor_by_gene$pearson_log2_tpm_plus1))),
  paste0("Median Pearson log2(TPM+1): ",
         median(cor_by_gene$pearson_log2_tpm_plus1, na.rm = TRUE)),
  paste0("Mean Pearson log2(TPM+1): ",
         mean(cor_by_gene$pearson_log2_tpm_plus1, na.rm = TRUE)),
  paste0("Median Spearman log2(TPM+1): ",
         median(cor_by_gene$spearman_log2_tpm_plus1, na.rm = TRUE)),
  paste0("Mean Spearman log2(TPM+1): ",
         mean(cor_by_gene$spearman_log2_tpm_plus1, na.rm = TRUE)),
  "",
  paste0("Saved overlap samples: ", sample_file),
  paste0("Saved merged values: ", merged_file),
  paste0("Saved per-gene correlations: ", cor_file),
  paste0("Saved boxplot: ", plot_file)
)
writeLines(summary_lines, summary_file)

message("Saved overlap samples: ", sample_file)
message("Saved merged values: ", merged_file)
message("Saved per-gene correlations: ", cor_file)
message("Saved boxplot: ", plot_file)
message("Saved summary: ", summary_file)

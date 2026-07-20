# Combine per-chromosome heart protein weights trained with:
#   moderate pruning 100kb r2=0.99, alpha=0.5, lambda=min.

library(data.table)
library(dplyr)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

results_dir <- file.path(paper_dir, "Results", "heart_protein_weights", "training_model_weights")
chr_list <- as.integer(strsplit(Sys.getenv("PWAS_CHR_LIST", unset = paste(1:22, collapse = ",")), ",")[[1]])
output_prefix <- Sys.getenv(
  "PWAS_OUTPUT_PREFIX",
  unset = "_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts"
)

weight_files <- file.path(
  results_dir,
  sprintf("weights_heart_protein_cardiomyocytes_other%s_chr%d.csv", output_prefix, chr_list)
)
skipped_files <- file.path(
  results_dir,
  sprintf("weights_heart_protein_cardiomyocytes_other_skipped%s_chr%d.csv", output_prefix, chr_list)
)
intercept_files <- file.path(
  results_dir,
  sprintf("intercepts_heart_protein_cardiomyocytes_other%s_chr%d.csv", output_prefix, chr_list)
)
model_term_files <- file.path(
  results_dir,
  sprintf("model_terms_heart_protein_cardiomyocytes_other%s_chr%d.csv", output_prefix, chr_list)
)

existing_weight_files <- weight_files[file.exists(weight_files)]
existing_skipped_files <- skipped_files[file.exists(skipped_files)]
existing_intercept_files <- intercept_files[file.exists(intercept_files)]
existing_model_term_files <- model_term_files[file.exists(model_term_files)]

if (!length(existing_weight_files)) {
  stop("No per-chromosome moderate-pruned weight files found.", call. = FALSE)
}

weights <- bind_rows(lapply(existing_weight_files, fread))
intercepts <- bind_rows(lapply(existing_intercept_files, fread))
model_terms <- bind_rows(lapply(existing_model_term_files, fread))
skipped <- bind_rows(lapply(existing_skipped_files, function(path) {
  if (file.info(path)$size == 0) {
    return(data.frame())
  }
  fread(path)
}))

weights_out <- file.path(
  results_dir,
  paste0("weights_heart_protein_cardiomyocytes_other", output_prefix, ".csv")
)
skipped_out <- file.path(
  results_dir,
  paste0("weights_heart_protein_cardiomyocytes_other_skipped", output_prefix, ".csv")
)
intercepts_out <- file.path(
  results_dir,
  paste0("intercepts_heart_protein_cardiomyocytes_other", output_prefix, ".csv")
)
model_terms_out <- file.path(
  results_dir,
  paste0("model_terms_heart_protein_cardiomyocytes_other", output_prefix, ".csv")
)

fwrite(weights, weights_out)
fwrite(intercepts, intercepts_out)
fwrite(model_terms, model_terms_out)
fwrite(skipped, skipped_out)

snp_per_gene <- table(weights$gene_id)

cat("Combined chromosome weight files:", length(existing_weight_files), "\n")
cat("Combined chromosome intercept files:", length(existing_intercept_files), "\n")
cat("Combined chromosome model-term files:", length(existing_model_term_files), "\n")
cat("Rows:", nrow(weights), "\n")
cat("Genes:", length(unique(weights$gene_id)), "\n")
cat("Median SNP rows per gene:", median(snp_per_gene), "\n")
cat("Mean SNP rows per gene:", mean(snp_per_gene), "\n")
cat("Genes with 1 SNP:", sum(snp_per_gene == 1), "\n")
cat("Genes with >=5 SNP:", sum(snp_per_gene >= 5), "\n")
cat("Genes with >=10 SNP:", sum(snp_per_gene >= 10), "\n")
cat("Saved weights:", weights_out, "\n")
cat("Saved intercepts:", intercepts_out, "\n")
cat("Saved model terms:", model_terms_out, "\n")
cat("Saved skipped:", skipped_out, "\n")

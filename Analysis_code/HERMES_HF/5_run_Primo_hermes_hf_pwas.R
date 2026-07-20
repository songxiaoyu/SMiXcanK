# Run Primo pattern inference on HERMES_HF HF PWAS S-MiXcan results.

library(data.table)
library(dplyr)
library(SMiXcan)

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hf_pwas",
  "hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
results_dir <- file.path(workspace_dir, "hermes_hf_result")
data_dir <- file.path(paper_dir, "Data")

result_tag <- Sys.getenv("HERMES_HF_RESULT_TAG", unset = "fixed_0p200")
result_suffix <- if (nzchar(result_tag)) paste0("_", result_tag) else ""
combined_path <- file.path(results_dir, paste0("hermes_hf_result_pwas", result_suffix, ".csv"))
p_join_col <- "p_join"
if (!file.exists(combined_path)) {
  stop("Missing HERMES_HF HF PWAS result file: ", combined_path, call. = FALSE)
}

combined <- read.csv(combined_path)
if (!p_join_col %in% names(combined)) {
  stop("Missing requested HERMES_HF p-value column: ", p_join_col, call. = FALSE)
}
if (p_join_col != "p_join") {
  combined$p_join <- combined[[p_join_col]]
}

ensembl_ref <- read.csv(file.path(data_dir, "ensembl38.txt"))
setDT(combined)
combined[, gene_id_clean := sub("\\..*$", "", gene_id)]
setDT(ensembl_ref)

# Add cytoband/gene-symbol annotation for reporting. Keep gene_id as the unique
# model identifier because display gene symbols can be ambiguous in protein
# weight analyses.
ensembl_cyto <- unique(
  ensembl_ref[, .(
    ENSG = sub("\\..*$", "", Gene.stable.ID),
    CYTOBAND = Karyotype.band,
    gene_name_ref = Gene.name
  )],
  by = "ENSG"
)

combined <- ensembl_cyto[combined, on = .(ENSG = gene_id_clean)]
combined[, ENSG := NULL]
combined[, gene_name_final := fifelse(
  !is.na(gene_name) & gene_name != "",
  as.character(gene_name),
  as.character(gene_name_ref)
)]
setcolorder(
  combined,
  c(
    "gene_name_final", "gene_id", "CYTOBAND",
    setdiff(names(combined), c("gene_name_final", "gene_id", "CYTOBAND"))
  )
)

combined <- as.data.frame(combined)
primo <- infer_celltype_patterns(
  merged = combined,
  pvals_names = c("p_cardiomyocytes", "p_other"),
  p_join_name = "p_join",
  type_col = "type",
  gene_id_col = "gene_id"
)

out <- primo$out
# The annotated result should expose one display gene_name column only. Helper
# columns used during annotation are dropped to avoid duplicated gene-name
# columns in downstream result tables.
out_clean <- out %>%
  select(-any_of(c("gene_name", "gene_name_ref", "gene_id_clean"))) %>%
  rename(gene_name = gene_name_final)

annotated_path <- file.path(results_dir, paste0("hermes_hf_result_pwas", result_suffix, "_annotated.csv"))
write.csv(out_clean, annotated_path, row.names = FALSE)

out_rename <- out_clean %>%
  rename(
    model = type,
    Z_cardiomyocyte = Z_cardiomyocytes,
    p_cardiomyocyte = p_cardiomyocytes,
    p_joint = p_join,
    fdr_p_joint = fdr_p_join,
    prob_00 = post_00,
    prob_10 = post_10,
    prob_01 = post_01,
    prob_11 = post_11
  ) %>%
  select(
    gene_name, gene_id, chr, CYTOBAND, model, input_snp_num,
    Z_cardiomyocyte, p_cardiomyocyte, Z_other, p_other,
    p_joint, fdr_p_joint,
    prob_00, prob_10, prob_01, prob_11,
    MAP_pattern_nonnull
  )

table_path <- file.path(results_dir, paste0("hermes_hf_table_pwas", result_suffix, ".csv"))
write.csv(out_rename, table_path, row.names = FALSE)

final_table_path <- file.path(workspace_dir, "FINAL_HERMES_HF_PWAS_TABLE.csv")
write.csv(out_rename, final_table_path, row.names = FALSE)

write_final_readme <- function(tab) {
  lines <- c(
    "# Final HERMES HF PWAS Report",
    "",
    "This workspace root contains the final deliverables for the active HERMES HF PWAS analysis.",
    "",
    "## Final files",
    "",
    "- `FINAL_HERMES_HF_PWAS_TABLE.csv`: clean final result table for reporting and plotting.",
    "- `FINAL_HERMES_HF_PWAS_QQ_COUNT_VENN.pdf`: QQ plot plus count-based Venn figure.",
    "- `FINAL_HERMES_HF_PWAS_QQ_GENE_VENN.pdf`: QQ plot plus gene-name Venn figure.",
    "",
    "## Analysis settings",
    "",
    "- GWAS: HERMES HF",
    "- Heart protein weights: moderate 100kb, r2 0.99, alpha 0.5, lambda.min",
    "- Association regularization: fixed scale 0.2",
    "",
    "## Result summary",
    "",
    paste0("- Total tested genes/proteins: ", nrow(tab)),
    paste0("- FDR < 0.05: ", sum(tab$fdr_p_joint < 0.05, na.rm = TRUE)),
    paste0("- FDR < 0.10: ", sum(tab$fdr_p_joint < 0.10, na.rm = TRUE)),
    paste0("- Minimum joint p-value: ", format(min(tab$p_joint, na.rm = TRUE), scientific = TRUE)),
    "",
    "Original pipeline outputs are still preserved in `hermes_hf_result/`."
  )
  writeLines(lines, file.path(workspace_dir, "FINAL_REPORT_README.md"))
}
write_final_readme(out_rename)

cat("Input result:", combined_path, "\n")
cat("p_join column:", p_join_col, "\n")
cat("Annotated output:", annotated_path, "\n")
cat("Table output:", table_path, "\n")
cat("Final table:", final_table_path, "\n")

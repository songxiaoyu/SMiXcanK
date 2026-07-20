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

make_pvar_map <- function(pvar_file) {
  pvar <- fread(pvar_file, skip = "#CHROM", header = TRUE,
                select = c("#CHROM", "POS", "ID", "REF", "ALT"))
  setnames(pvar, "#CHROM", "CHR")
  pvar[, varID := paste(CHR, POS, REF, ALT, sep = ":")]
  pvar[, rev_varID := paste(CHR, POS, ALT, REF, sep = ":")]
  pvar
}

choose_dosage_col <- function(varID, rsid, dosed, ref, eff, raw_names) {
  direct <- paste0(varID, "_", dosed)
  if (direct %in% raw_names) return(direct)
  candidates <- unique(na.omit(c(
    paste0(rsid, "_", dosed),
    paste0(rsid, "_", ref),
    paste0(rsid, "_", eff)
  )))
  hit <- candidates[candidates %in% raw_names]
  if (length(hit)) hit[1] else NA_character_
}

predict_original110 <- function(
  chromosomes,
  dosage_dir,
  pvar_dir,
  smixcan_weights,
  predixcan_weights,
  intercepts_smixcan,
  intercepts_predixcan,
  pi_file,
  out_dir
) {
  pi <- fread(pi_file)
  pi_id_col <- setdiff(names(pi), c("Cardiomyocytes", "Other"))[1]
  setnames(pi, pi_id_col, "pi_sample_id")
  pi[, donor_id := donor_id(pi_sample_id)]
  pi <- pi[, .(
    pi_cardiomyocytes = mean(Cardiomyocytes, na.rm = TRUE),
    pi_other = mean(Other, na.rm = TRUE)
  ), by = donor_id]

  pred_files <- character()
  match_reports <- list()

  for (chr in chromosomes) {
    chr_value <- chr
    message("===== Predicting original 110 chr", chr, " =====")
    dosage_file <- file.path(
      dosage_dir,
      paste0("chr", chr, "_dosage_nomiss_LDpruned_100kb_1_r2_0.99.raw")
    )
    pvar_file <- file.path(pvar_dir, paste0("GTEx_EA_chr", chr, "_nomiss.pvar"))
    if (!file.exists(dosage_file)) stop("Missing dosage file: ", dosage_file)
    if (!file.exists(pvar_file)) stop("Missing pvar file: ", pvar_file)

    pvar <- make_pvar_map(pvar_file)
    raw_header <- names(fread(dosage_file, nrows = 0))
    id_cols <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"),
                         raw_header)

    sm_w <- smixcan_weights[as.integer(chr) == chr_value]
    px_w <- predixcan_weights[as.integer(chr) == chr_value]
    int_sm <- intercepts_smixcan[as.integer(chr) == chr_value]
    int_px <- intercepts_predixcan[as.integer(chr) == chr_value]

    sm_w <- merge(sm_w, pvar[, .(varID, ID)], by = "varID",
                  all.x = TRUE, sort = FALSE)
    missing_id <- is.na(sm_w$ID)
    if (any(missing_id)) {
      sm_w <- merge(sm_w, pvar[, .(rev_varID, ID_rev = ID)],
                    by.x = "varID", by.y = "rev_varID",
                    all.x = TRUE, sort = FALSE)
      sm_w[missing_id & !is.na(ID_rev), ID := ID_rev]
      sm_w[, ID_rev := NULL]
    }

    px_w <- merge(px_w, pvar[, .(varID, ID)], by = "varID",
                  all.x = TRUE, sort = FALSE)
    missing_id <- is.na(px_w$ID)
    if (any(missing_id)) {
      px_w <- merge(px_w, pvar[, .(rev_varID, ID_rev = ID)],
                    by.x = "varID", by.y = "rev_varID",
                    all.x = TRUE, sort = FALSE)
      px_w[missing_id & !is.na(ID_rev), ID := ID_rev]
      px_w[, ID_rev := NULL]
    }

    sm_w[, dosage_col := mapply(
      choose_dosage_col, varID, ID, dosed_allele, ref_allele, eff_allele,
      MoreArgs = list(raw_names = raw_header), USE.NAMES = FALSE
    )]
    px_w[, dosage_col := mapply(
      choose_dosage_col, varID, ID, dosed_allele, ref_allele, ref_allele,
      MoreArgs = list(raw_names = raw_header), USE.NAMES = FALSE
    )]

    match_reports[[as.character(chr)]] <- rbindlist(list(
      data.table(method = "SMiXcan", chr = chr_value,
                 n_weight_rows = nrow(sm_w),
                 n_with_dosage_col = sum(!is.na(sm_w$dosage_col)),
                 n_missing_dosage_col = sum(is.na(sm_w$dosage_col))),
      data.table(method = "PrediXcan", chr = chr_value,
                 n_weight_rows = nrow(px_w),
                 n_with_dosage_col = sum(!is.na(px_w$dosage_col)),
                 n_missing_dosage_col = sum(is.na(px_w$dosage_col)))
    ))

    use_cols <- unique(c(sm_w$dosage_col[!is.na(sm_w$dosage_col)],
                         px_w$dosage_col[!is.na(px_w$dosage_col)]))
    dosage <- fread(dosage_file, select = unique(c(id_cols, use_cols)))
    sample_dt <- dosage[, .(FID, IID)]
    sample_dt[, sample_id := IID]
    sample_dt[, donor_id := donor_id(sample_id)]
    n_samples <- nrow(sample_dt)

    dosage_mat <- as.matrix(dosage[, ..use_cols])
    storage.mode(dosage_mat) <- "numeric"
    colnames(dosage_mat) <- use_cols

    pred_list <- vector("list", nrow(int_sm))
    for (i in seq_len(nrow(int_sm))) {
      gene <- int_sm$gene_id[i]
      gene_sm <- sm_w[gene_id == gene & !is.na(dosage_col)]
      gene_px <- px_w[gene_id == gene & !is.na(dosage_col)]

      pred_cardio <- rep(int_sm$intercept_cardiomyocytes[i], n_samples)
      pred_other <- rep(int_sm$intercept_other[i], n_samples)
      pred_predixcan <- rep(
        int_px[gene_id == gene, intercept_tissue][1],
        n_samples
      )

      if (nrow(gene_sm)) {
        x <- dosage_mat[, gene_sm$dosage_col, drop = FALSE]
        pred_cardio <- pred_cardio + as.numeric(x %*% gene_sm$weight_cardiomyocytes)
        pred_other <- pred_other + as.numeric(x %*% gene_sm$weight_other)
      }
      if (nrow(gene_px)) {
        x <- dosage_mat[, gene_px$dosage_col, drop = FALSE]
        pred_predixcan <- pred_predixcan + as.numeric(x %*% gene_px$weight_tissue)
      }

      pred_dt <- data.table(
        donor_id = sample_dt$donor_id,
        sample_id = sample_dt$sample_id,
        gene_id = gene,
        gene_name = int_sm$gene_name[i],
        chr = chr_value,
        pred_cardiomyocytes = pred_cardio,
        pred_other = pred_other,
        pred_tissue_predixcan = pred_predixcan
      )
      pred_dt <- merge(pred_dt, pi, by = "donor_id", all.x = TRUE, sort = FALSE)
      pred_dt[, pred_tissue_smixcan :=
                pi_cardiomyocytes * pred_cardiomyocytes +
                pi_other * pred_other]
      pred_list[[i]] <- pred_dt
    }

    pred_chr <- rbindlist(pred_list, use.names = TRUE)
    pred_file <- file.path(out_dir,
                           paste0("original110_predictions_chr", chr, ".csv"))
    fwrite(pred_chr, pred_file)
    pred_files <- c(pred_files, pred_file)

    message("Rows written chr", chr, ": ", nrow(pred_chr))
    rm(dosage, dosage_mat, pred_chr, pred_list)
    gc()
  }

  prediction <- rbindlist(lapply(pred_files, fread), use.names = TRUE)
  match_report <- rbindlist(match_reports, use.names = TRUE)
  fwrite(prediction, file.path(out_dir, "original110_predictions_all_chr.csv"))
  fwrite(match_report, file.path(out_dir, "original110_snp_matching_report_by_chr.csv"))
  list(prediction = prediction, match_report = match_report)
}

make_truth_tables <- function(project_dir, protein_file, expression_file, target_genes) {
  load(protein_file)
  protein <- as.data.table(Imputed_protein)
  ann_cols <- c("gene_name", "gene_id", "chr", "start", "end")
  protein_sample_cols <- setdiff(names(protein), ann_cols)
  protein <- protein[gene_id %in% target_genes]
  protein_long <- melt(
    protein,
    id.vars = intersect(ann_cols, names(protein)),
    measure.vars = protein_sample_cols,
    variable.name = "protein_sample_id",
    value.name = "true_protein"
  )
  protein_long[, protein_sample_id := normalize_sample_id(protein_sample_id)]
  protein_long[, donor_id := donor_id(protein_sample_id)]
  protein_long[, true_protein := as.numeric(true_protein)]
  protein_truth <- protein_long[, .(
    true_protein = mean(true_protein, na.rm = TRUE),
    protein_sample_ids = paste(unique(protein_sample_id), collapse = ";")
  ), by = .(gene_id, donor_id)]

  expr <- fread(expression_file, skip = 2, check.names = FALSE)
  expr <- expr[Name %in% target_genes]
  sample_cols <- setdiff(names(expr), c("Name", "Description"))
  expr_long <- melt(
    expr,
    id.vars = c("Name", "Description"),
    measure.vars = sample_cols,
    variable.name = "expression_sample_id",
    value.name = "tpm"
  )
  setnames(expr_long, c("Name", "Description"),
           c("gene_id", "expression_gene_name"))
  expr_long[, donor_id := donor_id(expression_sample_id)]
  expr_long[, tpm := as.numeric(tpm)]
  expression_truth <- expr_long[, .(
    true_rna_tpm = mean(tpm, na.rm = TRUE),
    expression_sample_ids = paste(unique(expression_sample_id), collapse = ";")
  ), by = .(gene_id, donor_id)]
  expression_truth[, true_rna_log2_tpm_plus1 := log2(true_rna_tpm + 1)]

  list(protein = protein_truth, expression = expression_truth)
}

correlate_predictions <- function(prediction, protein_truth, expression_truth) {
  sm <- prediction[, .(
    donor_id, gene_id, gene_name, method = "SMiXcan",
    predicted = pred_tissue_smixcan
  )]
  px <- prediction[, .(
    donor_id, gene_id, gene_name, method = "PrediXcan",
    predicted = pred_tissue_predixcan
  )]
  pred_long <- rbindlist(list(sm, px), use.names = TRUE)

  protein_dat <- merge(pred_long, protein_truth,
                       by = c("gene_id", "donor_id"),
                       all = FALSE, sort = FALSE)
  protein_cor <- protein_dat[, .(
    truth = "true_protein",
    n = .N,
    pearson = safe_cor(predicted, true_protein, "pearson"),
    spearman = safe_cor(predicted, true_protein, "spearman")
  ), by = .(method, gene_id, gene_name)]

  expr_dat <- merge(pred_long, expression_truth,
                    by = c("gene_id", "donor_id"),
                    all = FALSE, sort = FALSE)
  expr_cor <- expr_dat[, .(
    truth = "true_rna_log2_tpm_plus1",
    n = .N,
    pearson = safe_cor(predicted, true_rna_log2_tpm_plus1, "pearson"),
    spearman = safe_cor(predicted, true_rna_log2_tpm_plus1, "spearman")
  ), by = .(method, gene_id, gene_name)]

  list(
    cor = rbindlist(list(protein_cor, expr_cor), use.names = TRUE),
    protein_dat = protein_dat,
    expr_dat = expr_dat
  )
}

project_dir <- get_env("PAPER_SMIXCAN_DIR",
                       "/Users/admin/Library/CloudStorage/Dropbox/Paper_SMiXcan")
out_dir <- get_env(
  "PWAS_ORIGINAL110_VALIDATION_OUT_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/original110_prediction_validation")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dosage_dir <- get_env(
  "PWAS_ORIGINAL110_DOSAGE_DIR",
  file.path(project_dir,
            "New generated files/codes/moderate_pruned_by_chr_100kb_1_r2_0.99")
)
pvar_dir <- get_env(
  "PWAS_ORIGINAL110_PVAR_DIR",
  file.path(project_dir, "New generated files/codes/by_chr_nomiss")
)
smix_prefix <- "_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts"
predix_prefix <- "_moderate_100kb_r2_0.99_alpha0.5_lambdamin_predixcan_with_intercepts"
weights_dir <- file.path(project_dir,
                         "Results/heart_protein_weights/training_model_weights")
smix_weights_file <- file.path(
  weights_dir,
  paste0("weights_heart_protein_cardiomyocytes_other", smix_prefix, ".csv")
)
smix_intercepts_file <- file.path(
  weights_dir,
  paste0("intercepts_heart_protein_cardiomyocytes_other", smix_prefix, ".csv")
)
predix_weights_file <- file.path(
  weights_dir,
  paste0("predixcan_tissue_weights_heart_protein", predix_prefix, ".csv")
)
predix_intercepts_file <- file.path(
  weights_dir,
  paste0("intercepts_heart_protein_cardiomyocytes_other", predix_prefix, ".csv")
)
pi_file <- file.path(project_dir, "Heart/GTEx_Pi_Estimate/RNA_all_2celltypes_pi.tsv")
protein_file <- file.path(
  project_dir,
  "Heart/GTEx_Pi_Estimate/Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData"
)
expression_file <- file.path(project_dir,
                             "Heart/Data/gene_tpm_2022-06-06_v10_heart_left_ventricle.gct")
chromosomes <- as.integer(strsplit(get_env("PWAS_CHR", paste(1:22, collapse = ",")),
                                   ",", fixed = TRUE)[[1]])
chromosomes <- chromosomes[!is.na(chromosomes)]

message("Output dir: ", out_dir)
message("Original 110 dosage dir: ", dosage_dir)

load(protein_file)
protein_gene_ids <- unique(as.data.table(Imputed_protein)$gene_id)
rm(Imputed_protein)
expr_gene_ids <- fread(expression_file, skip = 2, select = "Name")$Name
truth_overlap_genes <- intersect(protein_gene_ids, expr_gene_ids)
message("Genes with both true protein and true RNA available: ",
        length(truth_overlap_genes))

smix_weights <- fread(smix_weights_file)
predix_weights <- fread(predix_weights_file)
intercepts_smixcan <- fread(smix_intercepts_file)
intercepts_predixcan <- fread(predix_intercepts_file)

smix_weights <- smix_weights[gene_id %in% truth_overlap_genes]
predix_weights <- predix_weights[gene_id %in% truth_overlap_genes]
intercepts_smixcan <- intercepts_smixcan[gene_id %in% truth_overlap_genes]
intercepts_predixcan <- intercepts_predixcan[gene_id %in% truth_overlap_genes]
message("Validation genes retained in SMiXcan intercepts: ",
        uniqueN(intercepts_smixcan$gene_id))
message("Validation genes retained in PrediXcan intercepts: ",
        uniqueN(intercepts_predixcan$gene_id))

pred_result <- predict_original110(
  chromosomes = chromosomes,
  dosage_dir = dosage_dir,
  pvar_dir = pvar_dir,
  smixcan_weights = smix_weights,
  predixcan_weights = predix_weights,
  intercepts_smixcan = intercepts_smixcan,
  intercepts_predixcan = intercepts_predixcan,
  pi_file = pi_file,
  out_dir = out_dir
)

truth <- make_truth_tables(
  project_dir = project_dir,
  protein_file = protein_file,
  expression_file = expression_file,
  target_genes = unique(pred_result$prediction$gene_id)
)
cors <- correlate_predictions(pred_result$prediction, truth$protein, truth$expression)

cor_file <- file.path(out_dir, "original110_prediction_truth_correlation_by_gene.csv")
protein_merged_file <- file.path(out_dir, "original110_prediction_true_protein_merged.csv")
rna_merged_file <- file.path(out_dir, "original110_prediction_true_rna_merged.csv")
summary_file <- file.path(out_dir, "original110_prediction_validation_summary.txt")
plot_file <- file.path(out_dir, "original110_prediction_truth_correlation_boxplot.png")

fwrite(cors$cor, cor_file)
fwrite(cors$protein_dat, protein_merged_file)
fwrite(cors$expr_dat, rna_merged_file)

png(plot_file, width = 1800, height = 1400, res = 180)
boxplot(
  pearson ~ interaction(method, truth, sep = "\n"),
  data = cors$cor[!is.na(pearson)],
  col = c("#4C78A8", "#F58518", "#4C78A8", "#F58518"),
  ylab = "Pearson correlation",
  xlab = "",
  main = "Original 110: Predicted Protein vs Truth"
)
abline(h = 0, lty = 2, col = "gray40")
dev.off()

summary_dt <- cors$cor[, .(
  genes = .N,
  genes_nonmissing = sum(!is.na(pearson)),
  median_pearson = median(pearson, na.rm = TRUE),
  mean_pearson = mean(pearson, na.rm = TRUE),
  median_spearman = median(spearman, na.rm = TRUE),
  mean_spearman = mean(spearman, na.rm = TRUE)
), by = .(method, truth)]

summary_lines <- c(
  paste0("Original 110 validation date: ", Sys.time()),
  paste0("Dosage dir: ", dosage_dir),
  paste0("SMiXcan weights: ", smix_weights_file),
  paste0("PrediXcan weights: ", predix_weights_file),
  paste0("Prediction rows: ", nrow(pred_result$prediction)),
  paste0("Predicted donors: ", uniqueN(pred_result$prediction$donor_id)),
  paste0("Predicted genes: ", uniqueN(pred_result$prediction$gene_id)),
  paste0("Prediction/true protein merged rows: ", nrow(cors$protein_dat)),
  paste0("Prediction/true protein donors: ", uniqueN(cors$protein_dat$donor_id)),
  paste0("Prediction/true RNA merged rows: ", nrow(cors$expr_dat)),
  paste0("Prediction/true RNA donors: ", uniqueN(cors$expr_dat$donor_id)),
  "",
  "SNP matching:",
  capture.output(print(pred_result$match_report[, .(
    n_weight_rows = sum(n_weight_rows),
    n_with_dosage_col = sum(n_with_dosage_col),
    n_missing_dosage_col = sum(n_missing_dosage_col)
  ), by = method])),
  "",
  "Correlation summary:",
  capture.output(print(summary_dt)),
  "",
  paste0("Saved correlations: ", cor_file),
  paste0("Saved true protein merged data: ", protein_merged_file),
  paste0("Saved true RNA merged data: ", rna_merged_file),
  paste0("Saved boxplot: ", plot_file)
)
writeLines(summary_lines, summary_file)

message("Saved correlations: ", cor_file)
message("Saved summary: ", summary_file)

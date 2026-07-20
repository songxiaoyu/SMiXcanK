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

dosage_dir <- get_env(
  "PWAS_DOSAGE_DIR",
  file.path(project_dir, "New generated files/codes/by_chr_heart_left_ventricle_351_samples_variant_id")
)

weights_file <- get_env(
  "PWAS_WEIGHTS_FILE",
  file.path(project_dir,
            "Results/heart_protein_weights/training_model_weights",
            "weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts.csv")
)

intercepts_file <- get_env(
  "PWAS_INTERCEPTS_FILE",
  file.path(project_dir,
            "Results/heart_protein_weights/training_model_weights",
            "intercepts_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts.csv")
)

snp_map_file <- get_env(
  "PWAS_SNP_MAP_FILE",
  file.path(project_dir, "Data/plink_snplist_by_gene/ALL_genes_snplist_with_rsids.csv")
)

out_dir <- get_env(
  "PWAS_PRED_OUT_DIR",
  file.path(project_dir,
            "Results/heart_protein_weights/predicted_protein_HLV351_with_intercepts")
)

chromosomes <- as.integer(strsplit(get_env("PWAS_CHR", paste(1:22, collapse = ",")),
                                   ",", fixed = TRUE)[[1]])
chromosomes <- chromosomes[!is.na(chromosomes)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Project dir: ", project_dir)
message("Dosage dir: ", dosage_dir)
message("Weights: ", weights_file)
message("Intercepts: ", intercepts_file)
message("SNP map: ", snp_map_file)
message("Output dir: ", out_dir)
message("Chromosomes: ", paste(chromosomes, collapse = ","))

weights <- fread(weights_file)
intercepts <- fread(intercepts_file)
snp_map <- fread(snp_map_file, select = c("key", "rsid"))
snp_map <- unique(snp_map[!is.na(key) & !is.na(rsid) & nzchar(key) & nzchar(rsid)],
                  by = "key")

required_weight_cols <- c(
  "gene_id", "gene_name", "varID", "chr", "ref_allele", "eff_allele",
  "dosed_allele", "weight_cardiomyocytes", "weight_other"
)
missing_weight_cols <- setdiff(required_weight_cols, names(weights))
if (length(missing_weight_cols)) {
  stop("Weights file is missing columns: ", paste(missing_weight_cols, collapse = ", "))
}

required_intercept_cols <- c(
  "gene_id", "gene_name", "chr", "intercept_cardiomyocytes", "intercept_other"
)
missing_intercept_cols <- setdiff(required_intercept_cols, names(intercepts))
if (length(missing_intercept_cols)) {
  stop("Intercepts file is missing columns: ", paste(missing_intercept_cols, collapse = ", "))
}

reverse_variant_key <- function(x) {
  parts <- tstrsplit(x, ":", fixed = TRUE)
  if (length(parts) != 4L) {
    return(rep(NA_character_, length(x)))
  }
  paste(parts[[1]], parts[[2]], parts[[4]], parts[[3]], sep = ":")
}

pick_dosage_column <- function(rsid, dosed, eff, ref, raw_names) {
  candidates <- unique(na.omit(c(
    paste0(rsid, "_", dosed),
    paste0(rsid, "_", eff),
    paste0(rsid, "_", ref)
  )))
  hit <- candidates[candidates %in% raw_names]
  if (length(hit)) hit[1] else NA_character_
}

prediction_files <- character()
match_reports <- list()

for (chr in chromosomes) {
  chr_value <- chr
  message("===== Predicting chr", chr, " =====")
  dosage_file <- file.path(dosage_dir, paste0("chr", chr, "_HLV351_dosage_nomiss.raw"))
  if (!file.exists(dosage_file)) {
    dosage_file_variant_id <- file.path(
      dosage_dir,
      paste0("chr", chr, "_HLV351_dosage_nomiss_variant_id.raw")
    )
    if (file.exists(dosage_file_variant_id)) {
      dosage_file <- dosage_file_variant_id
    }
  }
  if (!file.exists(dosage_file)) {
    stop("Missing dosage file: ", dosage_file)
  }

  intercept_chr <- intercepts[as.integer(chr) == chr_value]
  if (!nrow(intercept_chr)) {
    warning("No intercept rows for chr", chr, "; skipping")
    next
  }

  weights_chr <- weights[as.integer(chr) == chr_value]
  weights_chr[, reverse_varID := reverse_variant_key(varID)]

  weights_chr <- merge(weights_chr, snp_map, by.x = "varID", by.y = "key",
                       all.x = TRUE, sort = FALSE)
  missing_rsid <- is.na(weights_chr$rsid)
  if (any(missing_rsid)) {
    rev_map <- copy(snp_map)
    setnames(rev_map, c("key", "rsid"), c("reverse_varID", "rsid_reverse"))
    weights_chr <- merge(weights_chr, rev_map, by = "reverse_varID",
                         all.x = TRUE, sort = FALSE)
    weights_chr[missing_rsid & !is.na(rsid_reverse), rsid := rsid_reverse]
    weights_chr[, rsid_reverse := NULL]
  }

  raw_header <- names(fread(dosage_file, nrows = 0))
  id_cols <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), raw_header)

  if (nrow(weights_chr)) {
    weights_chr[, direct_dosage_col := paste0(varID, "_", dosed_allele)]
    weights_chr[, dosage_col := fifelse(direct_dosage_col %in% raw_header,
                                         direct_dosage_col,
                                         NA_character_)]
    needs_rsid_match <- is.na(weights_chr$dosage_col)
    if (any(needs_rsid_match)) {
      weights_chr[needs_rsid_match, dosage_col := mapply(
        pick_dosage_column,
        rsid = rsid,
        dosed = dosed_allele,
        eff = eff_allele,
        ref = ref_allele,
        MoreArgs = list(raw_names = raw_header),
        USE.NAMES = FALSE
      )]
    }
    weights_chr[, direct_dosage_col := NULL]
  } else {
    weights_chr[, dosage_col := character()]
  }

  match_report_chr <- weights_chr[, .(
    chr = chr_value,
    n_weight_rows = .N,
    n_with_rsid = sum(!is.na(rsid)),
    n_with_dosage_col = sum(!is.na(dosage_col)),
    n_missing_dosage_col = sum(is.na(dosage_col))
  )]
  match_reports[[as.character(chr)]] <- match_report_chr

  missing_rows <- weights_chr[is.na(dosage_col),
                              .(gene_id, gene_name, varID, chr, ref_allele,
                                eff_allele, dosed_allele, rsid)]
  if (nrow(missing_rows)) {
    fwrite(missing_rows,
           file.path(out_dir, paste0("missing_weight_snps_chr", chr, ".csv")))
  }

  use_weights <- weights_chr[!is.na(dosage_col)]
  dosage_cols <- unique(use_weights$dosage_col)
  select_cols <- unique(c(id_cols, dosage_cols))
  dosage <- fread(dosage_file, select = select_cols)

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

    pred_cardio <- rep(intercept_chr$intercept_cardiomyocytes[i], n_samples)
    pred_other <- rep(intercept_chr$intercept_other[i], n_samples)

    if (nrow(gene_weights)) {
      gene_cols <- gene_weights$dosage_col
      x <- dosage_mat[, gene_cols, drop = FALSE]
      pred_cardio <- pred_cardio +
        as.numeric(x %*% gene_weights$weight_cardiomyocytes)
      pred_other <- pred_other +
        as.numeric(x %*% gene_weights$weight_other)
    }

    pred_list[[i]] <- data.table(
      FID = sample_dt$FID,
      IID = sample_dt$IID,
      sample_id = sample_dt$sample_id,
      gene_id = gene,
      gene_name = intercept_chr$gene_name[i],
      chr = chr,
      pred_cardiomyocytes = pred_cardio,
      pred_other = pred_other
    )
  }

  pred_chr <- rbindlist(pred_list, use.names = TRUE)
  pred_file <- file.path(out_dir, paste0("predicted_heart_protein_HLV351_chr",
                                         chr, "_with_intercepts.csv"))
  fwrite(pred_chr, pred_file)
  prediction_files <- c(prediction_files, pred_file)

  message("Rows written chr", chr, ": ", nrow(pred_chr))
  message("Matched weight rows chr", chr, ": ",
          match_report_chr$n_with_dosage_col, " / ", match_report_chr$n_weight_rows)

  rm(dosage, dosage_mat, pred_chr, pred_list)
  gc()
}

match_report <- rbindlist(match_reports, use.names = TRUE, fill = TRUE)
fwrite(match_report, file.path(out_dir, "snp_matching_report_by_chr.csv"))

message("Combining chromosome predictions...")
pred_all <- rbindlist(lapply(prediction_files, fread), use.names = TRUE)
combined_file <- file.path(out_dir,
                           "predicted_heart_protein_HLV351_all_chr_with_intercepts.csv")
fwrite(pred_all, combined_file)

summary_file <- file.path(out_dir, "prediction_summary.txt")
summary_lines <- c(
  paste0("Prediction date: ", Sys.time()),
  paste0("Dosage dir: ", dosage_dir),
  paste0("Weights file: ", weights_file),
  paste0("Intercepts file: ", intercepts_file),
  paste0("SNP map file: ", snp_map_file),
  paste0("Combined prediction file: ", combined_file),
  paste0("Samples: ", uniqueN(pred_all$sample_id)),
  paste0("Genes: ", uniqueN(pred_all$gene_id)),
  paste0("Rows: ", nrow(pred_all)),
  paste0("Total weight rows: ", nrow(weights)),
  paste0("Matched weight rows: ", sum(match_report$n_with_dosage_col)),
  paste0("Missing dosage-column weight rows: ", sum(match_report$n_missing_dosage_col))
)
writeLines(summary_lines, summary_file)

message("Saved combined predictions: ", combined_file)
message("Saved matching report: ",
        file.path(out_dir, "snp_matching_report_by_chr.csv"))
message("Saved summary: ", summary_file)

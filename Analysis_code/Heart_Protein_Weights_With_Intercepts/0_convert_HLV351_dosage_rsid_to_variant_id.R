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
in_dir <- get_env(
  "PWAS_DOSAGE_IN_DIR",
  file.path(project_dir, "New generated files/codes/by_chr_heart_left_ventricle_351_samples")
)
out_dir <- get_env(
  "PWAS_DOSAGE_OUT_DIR",
  file.path(project_dir,
            "New generated files/codes/by_chr_heart_left_ventricle_351_samples_variant_id")
)
pvar_dir <- get_env(
  "PWAS_PVAR_DIR",
  file.path(project_dir, "New generated files/codes/by_chr_nomiss")
)
chromosomes <- as.integer(strsplit(get_env("PWAS_CHR", paste(1:22, collapse = ",")),
                                   ",", fixed = TRUE)[[1]])
chromosomes <- chromosomes[!is.na(chromosomes)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_mapping <- function(pvar_file) {
  pvar <- fread(pvar_file, skip = "#CHROM", header = TRUE,
                select = c("#CHROM", "POS", "ID", "REF", "ALT"))
  setnames(pvar, "#CHROM", "CHR")
  pvar <- pvar[!is.na(ID) & nzchar(ID)]
  pvar[, variant_id := paste(CHR, POS, REF, ALT, sep = ":")]

  ref_map <- pvar[, .(
    old_col = paste0(ID, "_", REF),
    new_col = paste0(variant_id, "_", REF)
  )]
  alt_map <- pvar[, .(
    old_col = paste0(ID, "_", ALT),
    new_col = paste0(variant_id, "_", ALT)
  )]
  unique(rbindlist(list(ref_map, alt_map)), by = "old_col")
}

for (chr in chromosomes) {
  message("===== Converting chr", chr, " =====")
  in_file <- file.path(in_dir, paste0("chr", chr, "_HLV351_dosage_nomiss.raw"))
  out_file <- file.path(out_dir, paste0("chr", chr, "_HLV351_dosage_nomiss_variant_id.raw"))
  report_file <- file.path(out_dir, paste0("chr", chr, "_HLV351_dosage_column_rename_report.csv"))
  pvar_file <- file.path(pvar_dir, paste0("GTEx_EA_chr", chr, "_nomiss.pvar"))

  if (!file.exists(in_file)) stop("Missing input dosage file: ", in_file)
  if (!file.exists(pvar_file)) stop("Missing pvar file: ", pvar_file)

  mapping <- make_mapping(pvar_file)
  dosage <- fread(in_file)
  old_names <- names(dosage)
  id_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  variant_cols <- setdiff(old_names, id_cols)

  rename_dt <- data.table(old_col = variant_cols)
  rename_dt <- merge(rename_dt, mapping, by = "old_col", all.x = TRUE, sort = FALSE)
  rename_dt[, renamed := !is.na(new_col)]
  rename_dt[is.na(new_col), new_col := old_col]

  setnames(dosage, old = rename_dt$old_col, new = rename_dt$new_col,
           skip_absent = TRUE)
  fwrite(dosage, out_file, sep = "\t")
  fwrite(rename_dt, report_file)

  message("Input variant columns: ", length(variant_cols))
  message("Renamed columns: ", sum(rename_dt$renamed))
  message("Unmapped columns kept unchanged: ", sum(!rename_dt$renamed))
  message("Saved: ", out_file)
  message("Saved report: ", report_file)
}

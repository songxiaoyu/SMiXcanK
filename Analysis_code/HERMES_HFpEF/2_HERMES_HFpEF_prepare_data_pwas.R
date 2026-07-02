# Prepare HERMES_HFpEF HFpEF GWAS for PWAS/S-MiXcan analysis.
#
# Inputs:
#   1. Heart protein weights from Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R
#   2. HERMES HFpEF GWAS lifted from hg19 to hg38 by 1_liftover_hermes_hfpef_gwas.py
#
# Output:
#   1. Per-chromosome SNP ID lists for 1000Genome LD extraction
#   2. Per-chromosome merged weight + GWAS RDS files for S-MiXcan

library(data.table)
library(dplyr)

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
weights_dir <- file.path(paper_dir, "Results", "heart_protein_weights", "training_model_weights")
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hfpef_pwas",
  "hermes_hfpef_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
gwas_file <- file.path(
  paper_dir,
  "Heart", "HERMES_HFpEF",
  "FORMAT-METAL_Pheno4_EUR_hg38_rsid.tsv.gz"
)
weights_file <- file.path(
  weights_dir,
  "weights_heart_protein_cardiomyocytes_other_moderate_100kb_r2_0.99_alpha0.5_lambdamin.csv"
)
chr_list <- 1:22

input_dir <- file.path(workspace_dir, "hermes_hfpef_input")
filtered_id_dir <- file.path(workspace_dir, "hermes_hfpef_filtered_id")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtered_id_dir, recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  sub("^chr", "", as.character(x))
}

normalize_varid <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x)
  x <- sub("_b38$", "", x)
  gsub("_", ":", x)
}

# Step 1. Read trained heart-protein weights.
mw_input <- fread(weights_file)
required_weight_cols <- c(
  "gene_id", "gene_name", "varID", "chr", "pos",
  "ref_allele", "eff_allele", "weight_cardiomyocytes", "weight_other", "type"
)
missing_weight_cols <- setdiff(required_weight_cols, names(mw_input))
if (length(missing_weight_cols)) {
  stop("Weights file is missing columns: ", paste(missing_weight_cols, collapse = ", "))
}

mw_input <- mw_input %>%
  mutate(
    CHR = paste0("chr", normalize_chr(chr)),
    varID = normalize_varid(varID)
  )
setDT(mw_input)
setkey(mw_input, CHR, varID)

# Step 2. Read lifted HERMES_HFpEF GWAS once and standardize columns to the
# allele names expected by the association code.
gwas <- fread(gwas_file)
required_gwas_cols <- c("CHR38", "POS38", "A1", "A2", "A1_beta", "se", "pval")
missing_gwas_cols <- setdiff(required_gwas_cols, names(gwas))
if (length(missing_gwas_cols)) {
  stop("HERMES_HFpEF GWAS is missing columns: ", paste(missing_gwas_cols, collapse = ", "))
}
if (!"rsid" %in% names(gwas)) {
  gwas[, rsid := NA_character_]
}

# HERMES_HFpEF A1_beta is the beta for A1. Use A1 as Effect.Gwas and A2 as baseline.
gwas <- gwas %>%
  mutate(
    CHR = paste0("chr", normalize_chr(CHR38)),
    POS = as.integer(POS38),
    Baseline.Gwas = as.character(A2),
    Effect.Gwas = as.character(A1),
    beta.Gwas = as.numeric(A1_beta),
    SE.Gwas = as.numeric(se),
    pval.Gwas = as.numeric(pval),
    rsid = as.character(rsid)
  ) %>%
  filter(!is.na(POS), !is.na(beta.Gwas), !is.na(SE.Gwas))
setDT(gwas)

# Step 3. For each chromosome, match weight SNPs to GWAS SNPs in both allele
# orientations. Reverse matches have beta flipped so the GWAS effect allele
# agrees with the model effect allele.
for (chr in chr_list) {
  cat("Processing chromosome", chr, "\n")

  chr_label <- paste0("chr", chr)
  gwas_chr <- gwas[CHR == chr_label]
  if (nrow(gwas_chr) == 0) {
    cat("  --> No GWAS rows, skipping\n")
    next
  }

  gwas_fwd <- copy(gwas_chr)
  gwas_fwd[, varID := paste0(chr, ":", POS, ":", Baseline.Gwas, ":", Effect.Gwas)]
  gwas_fwd[, flip := FALSE]

  gwas_rev <- copy(gwas_chr)
  gwas_rev[, varID := paste0(chr, ":", POS, ":", Effect.Gwas, ":", Baseline.Gwas)]
  gwas_rev[, flip := TRUE]

  gwas_long <- rbindlist(list(gwas_fwd, gwas_rev), use.names = TRUE, fill = TRUE)
  mw_gwas_input <- mw_input[CHR == chr_label][gwas_long, on = .(varID), nomatch = 0L]

  if (nrow(mw_gwas_input) == 0) {
    cat("  --> No overlapping weighted GWAS SNPs\n")
  } else {
    mw_gwas_input[flip == TRUE, beta.Gwas := -beta.Gwas]
    mw_gwas_input[flip == TRUE, c("Baseline.Gwas", "Effect.Gwas") := list(Effect.Gwas, Baseline.Gwas)]
  }

  # Step 3 writes both inputs needed downstream:
  #   1. SNP IDs for extracting 1000 Genomes EUR reference genotypes.
  #   2. Full merged weights + GWAS rows for S-MiXcan association.
  write.table(
    data.frame(mw_gwas_input$varID),
    file = file.path(filtered_id_dir, sprintf("hermes_hfpef_filtered_chr%d_gwas_id_pwas.txt", chr)),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(
    mw_gwas_input,
    file = file.path(input_dir, sprintf("chr%d_mw_gwas_input_hermes_hfpef_pwas.rds", chr))
  )
}

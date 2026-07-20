# Train 2-component PWAS MiXcan models on GTEx heart proteomics -----------------
#
# Goal:
#   Build protein prediction weights for a downstream PWAS/S-MiXcan analysis.
#   Each protein is modeled separately using cis genotype predictors, covariates,
#   and two cell-fraction components:
#     1. Cardiomyocytes
#     2. Other = Fibroblasts + Endothelial + Smooth muscle/pericyte + Immune
#
# Workflow steps:
#   Step 1. Set paths and required input files.
#   Step 2. Define helper functions for sample IDs, protein/cell-fraction data,
#           annotation, genotype dosage parsing, and output formatting.
#   Step 3. Load protein, cell-fraction, covariate, and EA-sample inputs.
#   Step 4. Align all training inputs by GTEx donor ID.
#   Step 5. Train one PWAS MiXcan model per protein:
#           - load matching GTEx heart genotype dosage
#           - keep cis SNPs within +/- 1 Mb when using chromosome-wide dosage
#           - remove incomplete samples and low-dosage SNPs
#           - fit MiXcan_train_K with K = 2 cell components
#   Step 6. Save nonzero SNP weights and a skipped-protein log.
#
# This script does not read GWAS summary statistics and does not run liftover.

# ------------------------------------------------------------------------------
# Step 1. Set paths and required input files
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(tibble)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)

heart_dir <- file.path(paper_dir, "Heart")
protein_dir <- file.path(heart_dir, "GTEx_Pi_Estimate")
data_dir <- file.path(paper_dir, "Data")
results_dir <- file.path(paper_dir, "Results", "heart_protein_weights", "training_model_weights")

protein_file <- file.path(
  protein_dir,
  "Imputed_Bulkprotein_GTEx.Proteomics.pQTL_Input.Heart_20250215.protein_normalized.RData"
)
pi_file <- file.path(protein_dir, "BayesDeBulk_pi.tsv")
covariate_file <- file.path(paper_dir, "New generated files", "covariate_EA_with_age.txt")
ea_keep_file <- file.path(protein_dir, "GTEx_heart_EA_subject_ids.txt")
ensembl_file <- file.path(data_dir, "ensembl38.txt")
rsid_annotation_file <- Sys.getenv(
  "PWAS_RSID_ANNOT_FILE",
  unset = file.path(
    paper_dir,
    "New generated files", "codes", "pruned_by_chr",
    "pruned_variants_500kb_1_r2_0.8.annotation.tsv"
  )
)

# Expected input: GTEx heart WGS dosage files exported by PLINK with --export A.
# Default location:
#   Paper_SMiXcan/New generated files/codes/by_chr_nomiss
# These are unpruned dosage files where only SNPs with any missing genotype have
# been removed. Create them with:
#   New generated files/codes/by_chr_nomiss is the no-missing dosage source used
#   before moderate LD pruning.
#
# To reproduce the archived pruned experiment, set PWAS_GENO_RAW_DIR to:
#   Paper_SMiXcan/New generated files/codes/pruned_by_chr
geno_raw_dir <- Sys.getenv(
  "PWAS_GENO_RAW_DIR",
  unset = file.path(paper_dir, "New generated files", "codes", "by_chr_nomiss")
)
chr_filter <- Sys.getenv("PWAS_CHR_FILTER", unset = "")
output_suffix <- Sys.getenv("PWAS_OUTPUT_SUFFIX", unset = "_unpruned_nomiss_with_intercepts")
output_prefix <- Sys.getenv(
  "PWAS_OUTPUT_PREFIX",
  unset = "_moderate_100kb_r2_0.99_alpha0.5_lambdamin_with_intercepts"
)
combine_only <- Sys.getenv("PWAS_COMBINE_ONLY", unset = "0") == "1"
max_proteins <- as.integer(Sys.getenv("PWAS_MAX_PROTEINS", unset = "0"))
weight_eps <- as.numeric(Sys.getenv("PWAS_WEIGHT_EPS", unset = "1e-8"))
mixcan_alpha <- as.numeric(Sys.getenv("PWAS_MIXCAN_ALPHA", unset = "0.5"))
lambda_choice <- Sys.getenv("PWAS_LAMBDA_CHOICE", unset = "1se")
if (!lambda_choice %in% c("1se", "min")) {
  stop("PWAS_LAMBDA_CHOICE must be '1se' or 'min'.", call. = FALSE)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

combine_chromosome_outputs <- function(results_dir, output_prefix, chr_list) {
  weight_files <- file.path(
    results_dir,
    sprintf("weights_heart_protein_cardiomyocytes_other%s_chr%d.csv",
            output_prefix, chr_list)
  )
  skipped_files <- file.path(
    results_dir,
    sprintf("weights_heart_protein_cardiomyocytes_other_skipped%s_chr%d.csv",
            output_prefix, chr_list)
  )
  intercept_files <- file.path(
    results_dir,
    sprintf("intercepts_heart_protein_cardiomyocytes_other%s_chr%d.csv",
            output_prefix, chr_list)
  )
  tissue_weight_files <- file.path(
    results_dir,
    sprintf("predixcan_tissue_weights_heart_protein%s_chr%d.csv",
            output_prefix, chr_list)
  )
  model_term_files <- file.path(
    results_dir,
    sprintf("model_terms_heart_protein_cardiomyocytes_other%s_chr%d.csv",
            output_prefix, chr_list)
  )

  existing_weight_files <- weight_files[file.exists(weight_files)]
  existing_skipped_files <- skipped_files[file.exists(skipped_files)]
  existing_intercept_files <- intercept_files[file.exists(intercept_files)]
  existing_tissue_weight_files <- tissue_weight_files[file.exists(tissue_weight_files)]
  existing_model_term_files <- model_term_files[file.exists(model_term_files)]

  if (!length(existing_weight_files)) {
    stop("No per-chromosome weight files found for prefix: ", output_prefix,
         call. = FALSE)
  }
  if (length(existing_intercept_files) != length(existing_weight_files)) {
    stop("Weight/intercept chromosome file count mismatch.", call. = FALSE)
  }
  if (length(existing_tissue_weight_files) != length(existing_weight_files)) {
    stop("Weight/PrediXcan tissue-weight chromosome file count mismatch.",
         call. = FALSE)
  }
  if (length(existing_model_term_files) != length(existing_weight_files)) {
    stop("Weight/model-term chromosome file count mismatch.", call. = FALSE)
  }

  weights <- bind_rows(lapply(existing_weight_files, fread))
  intercepts <- bind_rows(lapply(existing_intercept_files, fread))
  tissue_weights <- bind_rows(lapply(existing_tissue_weight_files, fread))
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
    paste0("weights_heart_protein_cardiomyocytes_other_skipped",
           output_prefix, ".csv")
  )
  intercepts_out <- file.path(
    results_dir,
    paste0("intercepts_heart_protein_cardiomyocytes_other",
           output_prefix, ".csv")
  )
  tissue_weights_out <- file.path(
    results_dir,
    paste0("predixcan_tissue_weights_heart_protein", output_prefix, ".csv")
  )
  model_terms_out <- file.path(
    results_dir,
    paste0("model_terms_heart_protein_cardiomyocytes_other",
           output_prefix, ".csv")
  )

  fwrite(weights, weights_out)
  fwrite(intercepts, intercepts_out)
  fwrite(tissue_weights, tissue_weights_out)
  fwrite(model_terms, model_terms_out)
  fwrite(skipped, skipped_out)

  snp_per_gene <- table(weights$gene_id)

  cat("Combined chromosome weight files:", length(existing_weight_files), "\n")
  cat("Combined chromosome intercept files:", length(existing_intercept_files), "\n")
  cat("Combined chromosome PrediXcan tissue-weight files:",
      length(existing_tissue_weight_files), "\n")
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
  cat("Saved PrediXcan tissue weights:", tissue_weights_out, "\n")
  cat("Saved model terms:", model_terms_out, "\n")
  cat("Saved skipped:", skipped_out, "\n")
}

if (combine_only) {
  chr_list <- as.integer(strsplit(
    Sys.getenv("PWAS_CHR_LIST", unset = paste(1:22, collapse = ",")),
    ",",
    fixed = TRUE
  )[[1]])
  chr_list <- chr_list[!is.na(chr_list)]
  combine_chromosome_outputs(results_dir, output_prefix, chr_list)
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# Step 2. Define helper functions
# ------------------------------------------------------------------------------

# Convert GTEx sample IDs from R-safe column names back to dash-separated IDs.
# Example: GTEX.QESD.0526.SM.DEJHB -> GTEX-QESD-0526-SM-DEJHB.
normalize_sample_id <- function(x) {
  gsub("\\.", "-", x)
}

# Reduce a GTEx sample ID to the donor ID used for matching genotype, covariate,
# protein, and cell-fraction records.
# Example: GTEX-QESD-0526-SM-DEJHB -> GTEX-QESD.
donor_id <- function(x) {
  x <- normalize_sample_id(x)
  vapply(strsplit(x, "-"), function(parts) {
    paste(parts[seq_len(min(2L, length(parts)))], collapse = "-")
  }, character(1))
}

required_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
  path
}

# Covariates are stored as one row per covariate and one column per sample.
# MiXcan needs one row per sample, so transpose the table and add donor IDs.
read_covariates <- function(path) {
  cov_wide <- fread(required_file(path, "Covariate file"), check.names = FALSE)
  cov_name_col <- names(cov_wide)[1]

  cov <- as.data.frame(t(as.data.frame(cov_wide[, -1, with = FALSE])))
  colnames(cov) <- cov_wide[[cov_name_col]]
  cov$SampleID <- rownames(cov)
  cov$DonorID <- donor_id(cov$SampleID)
  rownames(cov) <- NULL

  value_cols <- setdiff(names(cov), c("SampleID", "DonorID"))
  cov[value_cols] <- lapply(cov[value_cols], function(x) as.numeric(as.character(x)))
  cov
}

# The imputed protein file contains the training phenotypes. The original
# non-imputed protein file keeps more proteins but has missing values; this
# script uses Imputed_protein because model fitting expects complete y values.
load_protein_data <- function(path) {
  e <- new.env()
  load(required_file(path, "Protein RData file"), envir = e)
  if (!exists("Imputed_protein", envir = e)) {
    stop("Expected object Imputed_protein in ", path, call. = FALSE)
  }
  get("Imputed_protein", envir = e)
}

# Load five BayesDeBulk fractions, then collapse them to the requested 2-cell
# model: Cardiomyocytes versus all other heart cell types combined.
# The pi table has no sample ID column, so rows are matched to protein sample
# columns by order, following the pi-estimation workflow.
read_pi_data <- function(path, protein_sample_cols) {
  pi <- fread(required_file(path, "BayesDeBulk pi file"))
  source_cell_type_cols <- c(
    "Cardiomyocytes",
    "Fibroblasts",
    "Endothelial",
    "Smooth_muscle_cells_Pericyte",
    "Immune_cells"
  )
  missing_cols <- setdiff(source_cell_type_cols, names(pi))
  if (length(missing_cols)) {
    stop(
      "BayesDeBulk pi file is missing cell-type columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(pi) != length(protein_sample_cols)) {
    stop(
      "BayesDeBulk pi rows (", nrow(pi),
      ") do not match protein sample columns (", length(protein_sample_cols), ").",
      call. = FALSE
    )
  }

  pi$SampleID <- normalize_sample_id(protein_sample_cols)
  pi$DonorID <- donor_id(pi$SampleID)
  pi[, (source_cell_type_cols) := lapply(.SD, as.numeric), .SDcols = source_cell_type_cols]
  pi[, Other := Fibroblasts + Endothelial + Smooth_muscle_cells_Pericyte + Immune_cells]
  pi[, c("SampleID", "DonorID", "Cardiomyocytes", "Other"), with = FALSE]
}

read_ea_donors <- function(path) {
  keep <- fread(required_file(path, "EA donor keep file"), header = FALSE)
  unique(as.character(keep[[1]]))
}

# Keep the protein annotation needed for cis-window SNP selection and output.
annotation_from_protein <- function(protein) {
  ann <- protein[, c("gene_name", "gene_id", "chr", "start", "end")]
  names(ann) <- c("gene_name", "gene_id", "chr", "start", "end")
  ann$gene_id_no_version <- sub("\\..*$", "", ann$gene_id)
  ann
}

# Fill missing chromosome/start/end fields from ensembl38.txt when possible.
# Most records should already have annotation in the protein file; this is a
# fallback to reduce skipped proteins.
fill_missing_annotation <- function(ann, ensembl_file) {
  if (!file.exists(ensembl_file)) {
    return(ann)
  }

  ens <- fread(ensembl_file)
  setnames(
    ens,
    old = names(ens),
    new = make_clean_names(names(ens))
  )
  needed_cols <- c("gene_stable_id", "chromosome_scaffold_name",
                   "transcript_start_bp", "transcript_end_bp")
  if (!all(needed_cols %in% names(ens))) {
    return(ann)
  }

  ens <- ens[, ..needed_cols] %>%
    group_by(gene_stable_id) %>%
    summarise(
      chromosome_scaffold_name = first(chromosome_scaffold_name),
      transcript_start_bp = min(transcript_start_bp, na.rm = TRUE),
      transcript_end_bp = max(transcript_end_bp, na.rm = TRUE),
      .groups = "drop"
    )
  ann <- ann %>%
    left_join(
      ens,
      by = c("gene_id_no_version" = "gene_stable_id")
    ) %>%
    mutate(
      chr = ifelse(is.na(chr), chromosome_scaffold_name, chr),
      start = ifelse(is.na(start), transcript_start_bp, start),
      end = ifelse(is.na(end), transcript_end_bp, end)
    ) %>%
    select(gene_name, gene_id, gene_id_no_version, chr, start, end)
  ann
}

make_clean_names <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# Genotype input is chromosome-wide PLINK raw dosage:
#   chr<chr>_dosage_nomiss.raw
# The script filters each chromosome file to cis SNPs using the protein gene
# coordinates below.
candidate_raw_files <- function(target_id, chr) {
  parent_raw_dir <- dirname(geno_raw_dir)
  primary_matches <- Sys.glob(file.path(geno_raw_dir, sprintf("chr%s_dosage*.raw", chr)))
  search_dirs <- unique(c(
    geno_raw_dir,
    file.path(geno_raw_dir, "pruned_by_chr"),
    file.path(parent_raw_dir, "pruned_by_chr")
  ))
  raw_names <- c(
    sprintf("chr%s_dosage_nomiss.raw", chr),
    sprintf("chr%s_dosage.raw", chr),
    sprintf("chr%s_dosage_nomiss_LDpruned_500kb_1_r2_0.8.raw", chr)
  )
  unique(c(primary_matches, unlist(lapply(search_dirs, file.path, raw_names), use.names = FALSE)))
}

# PLINK --export A may produce columns like:
#   1:11260293:G:A_G
# or, when variant IDs are rsIDs:
#   rs10904045_C
# The first format carries its own coordinates. The second format needs the
# matching chromosome pvar file to recover chr/pos/ref/alt annotation.
candidate_pvar_files <- function(chr) {
  parent_raw_dir <- dirname(geno_raw_dir)
  unique(c(
    file.path(geno_raw_dir, sprintf("GTEx_EA_chr%s_nomiss.pvar", chr)),
    file.path(geno_raw_dir, sprintf("chr%s_hg38.pvar", chr)),
    file.path(parent_raw_dir, "pruned_by_chr", sprintf("GTEx_EA_chr%s_nomiss.pvar", chr)),
    file.path(parent_raw_dir, "pruned_by_chr", sprintf("chr%s_hg38.pvar", chr))
  ))
}

read_pvar_annotation <- function(chr) {
  pvar_file <- candidate_pvar_files(chr)
  pvar_file <- pvar_file[file.exists(pvar_file)][1]
  if (is.na(pvar_file)) {
    return(NULL)
  }

  pvar <- fread(pvar_file, skip = "#CHROM", check.names = FALSE)
  names(pvar) <- sub("^#", "", names(pvar))
  needed <- c("CHROM", "POS", "ID", "REF", "ALT")
  if (!all(needed %in% names(pvar))) {
    return(NULL)
  }

  pvar[, c("CHROM", "POS", "ID", "REF", "ALT"), with = FALSE] %>%
    rename(
      chr = CHROM,
      pos = POS,
      rsid = ID,
      ref_allele = REF,
      alt_allele = ALT
    ) %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      rsid = as.character(rsid),
      ref_allele = as.character(ref_allele),
      alt_allele = as.character(alt_allele)
    )
}

parse_position_variant_cols <- function(raw_cols) {
  var_id <- sub("_[^_]+$", "", raw_cols)
  dosed_allele <- sub("^.*_", "", raw_cols)
  parts <- strsplit(var_id, ":", fixed = TRUE)

  data.frame(
    varID = var_id,
    raw_col = raw_cols,
    chr = vapply(parts, function(x) x[1], character(1)),
    pos = as.integer(vapply(parts, function(x) x[2], character(1))),
    ref_allele = vapply(parts, function(x) x[3], character(1)),
    eff_allele = vapply(parts, function(x) x[4], character(1)),
    dosed_allele = dosed_allele,
    stringsAsFactors = FALSE
  )
}

read_rsid_position_annotation <- function(chr) {
  if (!nzchar(rsid_annotation_file) || !file.exists(rsid_annotation_file)) {
    return(NULL)
  }

  rsid_pos <- fread(rsid_annotation_file, check.names = FALSE)
  names(rsid_pos) <- make_clean_names(names(rsid_pos))
  rsid_col <- intersect(names(rsid_pos), c("rsid", "rs_id", "id"))[1]
  chr_col <- intersect(names(rsid_pos), c("chr38", "chrom38", "chromosome38", "chr", "chromosome"))[1]
  pos_col <- intersect(names(rsid_pos), c("pos38", "position38", "bp38", "pos", "position", "bp"))[1]
  ref_col <- intersect(names(rsid_pos), c("ref", "ref_allele", "reference_allele"))[1]
  alt_col <- intersect(names(rsid_pos), c("alt", "alt_allele", "effect_allele"))[1]
  dosage_col <- intersect(names(rsid_pos), c("dosage_col", "raw_col", "dosage_column"))[1]
  if (is.na(rsid_col) || is.na(chr_col) || is.na(pos_col)) {
    stop(
      "PWAS_RSID_ANNOT_FILE must contain rsid plus chromosome and position columns.",
      call. = FALSE
    )
  }

  out <- rsid_pos %>%
    transmute(
      rsid = as.character(.data[[rsid_col]]),
      chr = sub("^chr", "", as.character(.data[[chr_col]])),
      pos = as.integer(.data[[pos_col]]),
      ref_allele = if (!is.na(ref_col)) as.character(.data[[ref_col]]) else NA_character_,
      alt_allele = if (!is.na(alt_col)) as.character(.data[[alt_col]]) else NA_character_,
      annotation_raw_col = if (!is.na(dosage_col)) as.character(.data[[dosage_col]]) else NA_character_
    ) %>%
    filter(chr == sub("^chr", "", as.character(chr))) %>%
    distinct(rsid, .keep_all = TRUE)
  out
}

parse_rsid_variant_cols <- function(raw_cols, chr) {
  pvar <- read_pvar_annotation(chr)

  raw_map <- data.frame(
    rsid = sub("_[^_]+$", "", raw_cols),
    raw_col = raw_cols,
    dosed_allele = sub("^.*_", "", raw_cols),
    stringsAsFactors = FALSE
  )

  # First try direct matching when the pvar ID column is also rsID.
  if (!is.null(pvar)) {
    snp_annot <- raw_map %>%
      inner_join(pvar, by = "rsid")
  } else {
    snp_annot <- data.frame()
  }

  # Most current pvar files use coordinate IDs instead. In that case, use the
  # pruned annotation file in pruned_by_chr. If it has ref/alt, no pvar join is
  # needed; otherwise recover ref/alt from pvar by chr/pos.
  if (nrow(snp_annot) == 0) {
    rsid_pos <- read_rsid_position_annotation(chr)
    if (!is.null(rsid_pos)) {
      snp_annot <- raw_map %>%
        inner_join(rsid_pos, by = "rsid")
      if (
        nrow(snp_annot) > 0 &&
        !is.null(pvar) &&
        any(is.na(snp_annot$ref_allele) | is.na(snp_annot$alt_allele))
      ) {
        snp_annot <- snp_annot %>%
          select(-ref_allele, -alt_allele) %>%
          inner_join(
            pvar %>% select(chr, pos, ref_allele, alt_allele),
            by = c("chr", "pos")
          )
      }
    }
  }

  if (nrow(snp_annot) == 0) {
    return(data.frame())
  }

  snp_annot %>%
    mutate(
      alt_choices = strsplit(alt_allele, ",", fixed = TRUE),
      eff_allele = vapply(seq_along(alt_choices), function(i) {
        choices <- alt_choices[[i]]
        if (dosed_allele[i] %in% choices) {
          dosed_allele[i]
        } else {
          choices[1]
        }
      }, character(1)),
      allele_match = dosed_allele == ref_allele | dosed_allele == eff_allele,
      varID = paste(chr, pos, ref_allele, eff_allele, sep = ":")
    ) %>%
    arrange(raw_col, desc(allele_match)) %>%
    distinct(raw_col, .keep_all = TRUE) %>%
    select(varID, raw_col, chr, pos, ref_allele, eff_allele, dosed_allele) %>%
    as.data.frame()
}

parse_raw_variant_cols <- function(raw_cols, chr) {
  var_id <- sub("_[^_]+$", "", raw_cols)
  if (all(grepl("^[^:]+:[0-9]+:[^:]+:[^:]+$", var_id))) {
    return(parse_position_variant_cols(raw_cols))
  }
  parse_rsid_variant_cols(raw_cols, chr)
}

genotype_chr_cache <- new.env(parent = emptyenv())

# Read one chromosome's raw dosage and annotation once, then reuse it for all
# proteins on that chromosome.
read_genotype_for_chr <- function(target_id, chr) {
  cache_key <- as.character(chr)
  if (exists(cache_key, envir = genotype_chr_cache, inherits = FALSE)) {
    return(get(cache_key, envir = genotype_chr_cache, inherits = FALSE))
  }

  files <- candidate_raw_files(target_id, chr)
  raw_file <- files[file.exists(files)][1]
  if (is.na(raw_file)) {
    assign(cache_key, NULL, envir = genotype_chr_cache)
    return(NULL)
  }

  cat("Reading genotype raw file:", raw_file, "\n")
  raw <- fread(raw_file, check.names = FALSE)
  sample_col <- if ("IID" %in% names(raw)) "IID" else names(raw)[2]
  raw$DonorID <- donor_id(as.character(raw[[sample_col]]))

  variant_cols <- setdiff(
    names(raw),
    c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "DonorID")
  )
  snp_annot <- parse_raw_variant_cols(variant_cols, chr)
  out <- list(raw_file = raw_file, raw = raw, snp_annot = snp_annot)
  assign(cache_key, out, envir = genotype_chr_cache)
  out
}

# Select genotype dosages for one protein's +/- 1 Mb cis window.
read_genotype_for_target <- function(target_id, chr, start, end, window_bp = 1e6) {
  geno_chr <- read_genotype_for_chr(target_id, chr)
  if (is.null(geno_chr)) {
    return(NULL)
  }

  snp_annot <- geno_chr$snp_annot
  cis_start <- max(0, as.integer(start) - window_bp)
  cis_end <- as.integer(end) + window_bp
  keep <- snp_annot$chr == as.character(chr) &
    snp_annot$pos >= cis_start &
    snp_annot$pos <= cis_end
  snp_annot <- snp_annot[keep, , drop = FALSE]

  if (nrow(snp_annot) == 0) {
    return(NULL)
  }

  x <- as.data.frame(geno_chr$raw[, c("DonorID", snp_annot$raw_col), with = FALSE])
  list(raw_file = geno_chr$raw_file, x = x, snp_annot = snp_annot)
}

# Turn MiXcan_train_K output into the flat weight table expected by downstream
# PWAS/GWAS preparation code.
combine_weights <- function(fit, snp_annot, target, cell_type_cols) {
  weights <- as.data.frame(fit$W)
  weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
  colnames(weights) <- weight_cols
  weights$varID <- rownames(fit$W)

  out <- snp_annot %>%
    select(varID, chr, pos, ref_allele, eff_allele, dosed_allele) %>%
    left_join(weights, by = "varID") %>%
    mutate(
      gene_id = target$gene_id,
      gene_name = target$gene_name,
      type = fit$type
    ) %>%
    select(
      gene_id, gene_name, varID, chr, pos, ref_allele, eff_allele,
      dosed_allele, all_of(weight_cols), type
    )
  out
}

combine_tissue_weights <- function(fit, snp_annot, target) {
  beta_all <- fit$beta.all.models
  missing_cols <- setdiff("Tissue", colnames(beta_all))
  if (length(missing_cols)) {
    stop("Missing Tissue column in beta.all.models.", call. = FALSE)
  }
  tissue_weights <- data.frame(
    varID = rownames(beta_all),
    weight_tissue = as.numeric(beta_all[, "Tissue"]),
    stringsAsFactors = FALSE
  )
  tissue_weights <- tissue_weights[
    tissue_weights$varID %in% snp_annot$varID,
    ,
    drop = FALSE
  ]

  out <- snp_annot %>%
    select(varID, chr, pos, ref_allele, eff_allele, dosed_allele) %>%
    left_join(tissue_weights, by = "varID") %>%
    mutate(
      gene_id = target$gene_id,
      gene_name = target$gene_name,
      type = "PrediXcan"
    ) %>%
    select(
      gene_id, gene_name, varID, chr, pos, ref_allele, eff_allele,
      dosed_allele, weight_tissue, type
    )
  out
}

extract_intercepts <- function(fit, target, cell_type_cols) {
  intercept_idx <- match("Intercept", rownames(fit$beta.all.models))
  if (is.na(intercept_idx)) {
    intercept_idx <- 1L
  }
  intercept_row <- fit$beta.all.models[intercept_idx, , drop = FALSE]
  cell_cols <- paste0("Cell", seq_along(cell_type_cols))
  missing_cols <- setdiff(c("Tissue", cell_cols), colnames(intercept_row))
  if (length(missing_cols)) {
    stop("Missing intercept columns in beta.all.models: ",
         paste(missing_cols, collapse = ", "))
  }
  out <- data.frame(
    gene_id = target$gene_id,
    gene_name = target$gene_name,
    chr = target$chr,
    type = fit$type,
    intercept_tissue = as.numeric(intercept_row[, "Tissue"]),
    stringsAsFactors = FALSE
  )
  for (k in seq_along(cell_type_cols)) {
    out[[paste0("intercept_", make_clean_names(cell_type_cols[k]))]] <-
      as.numeric(intercept_row[, cell_cols[k]])
  }
  out
}

intercepts_as_model_terms <- function(intercepts, cell_type_cols) {
  weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
  out <- data.frame(
    gene_id = intercepts$gene_id,
    gene_name = intercepts$gene_name,
    term_type = "intercept",
    varID = "Intercept",
    chr = intercepts$chr,
    pos = NA_integer_,
    ref_allele = NA_character_,
    eff_allele = NA_character_,
    dosed_allele = NA_character_,
    type = intercepts$type,
    stringsAsFactors = FALSE
  )
  for (k in seq_along(cell_type_cols)) {
    out[[weight_cols[k]]] <- intercepts[[paste0("intercept_", make_clean_names(cell_type_cols[k]))]]
  }
  out[, c("gene_id", "gene_name", "term_type", "varID", "chr", "pos",
          "ref_allele", "eff_allele", "dosed_allele", weight_cols, "type")]
}

# Local variant of SMiXcan::MiXcan_train_K that can use lambda.min instead of
# the package default lambda.1se. This is useful for parameter diagnostics when
# lambda.1se makes the trained PWAS weights too sparse.
MiXcan_train_K_select_lambda <- function(y, x, cov = NULL, pi_k, xNameMatrix = NULL,
                                         yName = NULL, foldid = NULL, alpha = 0.5,
                                         lambda_choice = "1se") {
  y <- as.matrix(y)
  x <- as.matrix(x)
  pi_k <- as.matrix(pi_k)
  n <- nrow(x)
  p <- ncol(x)
  K <- ncol(pi_k)
  if (nrow(pi_k) != n) {
    stop("pi_k must have the same number of rows as x (samples).")
  }
  if (is.null(xNameMatrix)) {
    if (!is.null(colnames(x))) {
      xNameMatrix <- data.frame(SNP = colnames(x), stringsAsFactors = FALSE)
    } else {
      xNameMatrix <- data.frame(SNP = paste0("SNP", seq_len(p)), stringsAsFactors = FALSE)
    }
  }
  if (is.null(foldid)) {
    set.seed(1L)
    foldid <- sample(1:10, n, replace = TRUE)
  }

  pick_lambda <- function(cv_fit) {
    if (lambda_choice == "min") cv_fit$lambda.min else cv_fit$lambda.1se
  }
  cv_mse_at_lambda <- function(cv_fit, lambda_value) {
    idx <- which.min(abs(log(cv_fit$lambda) - log(lambda_value)))
    cv_fit$cvm[idx]
  }
  cv_r2_at_lambda <- function(cv_fit, lambda_value, y_centered) {
    y_var <- mean(as.numeric(y_centered)^2)
    if (is.na(y_var) || y_var <= 0) {
      return(NA_real_)
    }
    1 - cv_mse_at_lambda(cv_fit, lambda_value) / y_var
  }

  y <- scale(y, center = TRUE, scale = FALSE)
  x <- scale(x, center = TRUE, scale = FALSE)
  if (is.null(cov)) {
    xcov <- x
    pcov <- 0L
  } else {
    cov <- as.matrix(cov)
    cov <- scale(cov, center = TRUE, scale = FALSE)
    pcov <- ncol(cov)
    xcov <- cbind(x, cov)
  }

  ft00 <- glmnet::cv.glmnet(x = xcov, y = y, family = "gaussian",
                            foldid = foldid, alpha = alpha)
  lambda_tissue <- pick_lambda(ft00)
  ft0 <- glmnet::glmnet(x = xcov, y = y, family = "gaussian",
                        lambda = lambda_tissue, alpha = alpha)
  est.tissue <- c(ft0$a0, as.numeric(ft0$beta))

  C <- stats::contr.helmert(K)
  c_mat <- pi_k %*% C
  Z_list <- vector("list", length = K - 1L)
  for (m in seq_len(K - 1L)) {
    Z_list[[m]] <- c_mat[, m] * x
  }
  Z_block <- do.call(cbind, Z_list)
  if (is.null(cov)) {
    XX <- cbind(c_mat, x, Z_block)
  } else {
    XX <- cbind(c_mat, x, Z_block, cov)
  }

  n_c <- K - 1L
  n_Z <- p * (K - 1L)
  penalty.factor <- c(rep(0, n_c), rep(1, p + n_Z + pcov))
  ft11 <- glmnet::cv.glmnet(x = XX, y = y, family = "gaussian",
                            foldid = foldid, alpha = alpha,
                            penalty.factor = penalty.factor)
  lambda_cell <- pick_lambda(ft11)
  ft <- glmnet::glmnet(x = XX, y = y, family = "gaussian",
                       lambda = lambda_cell, alpha = alpha,
                       penalty.factor = penalty.factor)
  est <- c(ft$a0, as.numeric(ft$beta))

  idx_c <- 1:n_c
  idx_barb <- (n_c + 1):(n_c + p)
  idx_Z <- (n_c + p + 1):(n_c + p + n_Z)
  idx_cov <- if (pcov > 0L) (length(est) - pcov + 1):length(est) else integer(0)

  alpha_vec <- est[idx_c]
  b_bar <- est[idx_barb]
  d_mat <- matrix(est[idx_Z], nrow = p, ncol = K - 1L, byrow = FALSE)
  a_cell <- as.numeric(ft$a0 + C %*% alpha_vec)
  B_mat <- matrix(NA_real_, nrow = p, ncol = K)
  for (k in seq_len(K)) {
    B_mat[, k] <- b_bar + d_mat %*% C[k, ]
  }
  colnames(B_mat) <- paste0("Cell", seq_len(K))
  rownames(B_mat) <- xNameMatrix[[1]]

  if (suppressWarnings(all(d_mat == 0))) {
    Type <- if (suppressWarnings(all(b_bar == 0))) "NoPredictor" else "NonSpecific"
  } else {
    Type <- "CellTypeSpecific"
  }

  beta.SNP.by.cell <- vector("list", length = K)
  for (k in seq_len(K)) {
    beta.SNP.by.cell[[k]] <- data.frame(weight = B_mat[, k], stringsAsFactors = FALSE)
    rownames(beta.SNP.by.cell[[k]]) <- rownames(B_mat)
  }
  names(beta.SNP.by.cell) <- paste0("Cell", seq_len(K))

  n_rows <- 1 + p + pcov
  beta_cell_mat <- matrix(NA_real_, nrow = n_rows, ncol = K)
  rownames(beta_cell_mat) <- c("Intercept", xNameMatrix[[1]],
                              if (pcov > 0) paste0("cov", seq_len(pcov)) else NULL)
  colnames(beta_cell_mat) <- paste0("Cell", seq_len(K))
  for (k in seq_len(K)) {
    beta_cell_mat[1, k] <- a_cell[k]
    beta_cell_mat[1 + (1:p), k] <- B_mat[, k]
    if (pcov > 0L) {
      beta_cell_mat[1 + p + (1:pcov), k] <- if (length(idx_cov) > 0L) est[idx_cov] else 0
    }
  }
  beta_all_models <- cbind(Tissue = est.tissue, beta_cell_mat)
  rownames(beta_all_models) <- rownames(beta_cell_mat)

  list(type = Type, beta.SNP.by.cell = beta.SNP.by.cell,
       beta.all.models = beta_all_models, W = B_mat,
       glmnet.cell = ft, glmnet.tissue = ft0, yName = yName,
       xNameMatrix = xNameMatrix,
       lambda_cell = lambda_cell,
       lambda_tissue = lambda_tissue,
       cv_mse_cell = cv_mse_at_lambda(ft11, lambda_cell),
       cv_mse_tissue = cv_mse_at_lambda(ft00, lambda_tissue),
       cv_r2_cell = cv_r2_at_lambda(ft11, lambda_cell, y),
       cv_r2_tissue = cv_r2_at_lambda(ft00, lambda_tissue, y))
}

# ------------------------------------------------------------------------------
# Step 3. Load protein, cell-fraction, covariate, and EA-sample inputs
# ------------------------------------------------------------------------------

# Protein columns are full sample IDs, while covariates/genotypes are matched by
# donor ID. sample_map keeps both forms so we can read y from the protein matrix
# and align all other inputs by DonorID.
protein <- load_protein_data(protein_file)
protein_sample_cols <- names(protein)[-(1:5)]
sample_map <- data.frame(
  RawSampleID = protein_sample_cols,
  SampleID = normalize_sample_id(protein_sample_cols),
  DonorID = donor_id(protein_sample_cols),
  stringsAsFactors = FALSE
) %>%
  distinct(DonorID, .keep_all = TRUE)

pi_data <- read_pi_data(pi_file, protein_sample_cols) %>%
  distinct(DonorID, .keep_all = TRUE)
cov_data <- read_covariates(covariate_file) %>%
  distinct(DonorID, .keep_all = TRUE)
ea_donors <- read_ea_donors(ea_keep_file)

# ------------------------------------------------------------------------------
# Step 4. Align all training inputs by GTEx donor ID
# ------------------------------------------------------------------------------

# Final training sample set: donors with protein abundance, cell fractions, and
# covariates available. The covariate file is already EA-filtered; the separate
# EA keep file is applied only when it has enough overlap with these samples.
sample_map <- sample_map %>%
  inner_join(pi_data, by = "DonorID", suffix = c("", "_pi")) %>%
  inner_join(cov_data, by = "DonorID", suffix = c("", "_cov"))

ea_overlap <- intersect(sample_map$DonorID, ea_donors)
if (length(ea_overlap) >= 20) {
  sample_map <- sample_map %>%
    filter(DonorID %in% ea_overlap)
} else {
  warning(
    "EA keep file overlaps only ", length(ea_overlap),
    " aligned protein/pi/covariate samples; using covariate-file samples instead."
  )
}

if (nrow(sample_map) < 20) {
  stop("Fewer than 20 aligned protein/pi/covariate EA samples.", call. = FALSE)
}

protein_ann <- annotation_from_protein(protein)
protein_ann <- fill_missing_annotation(protein_ann, ensembl_file)
protein_ann <- protein_ann %>%
  distinct(gene_id, .keep_all = TRUE)

# Optional test/debug filter. Example:
#   PWAS_CHR_FILTER=10 Rscript Analysis_code/Heart_Protein_Weights/2_train_heart_protein_weights.R
if (nzchar(chr_filter)) {
  keep_chr <- protein_ann$chr == chr_filter |
    protein_ann$chr == paste0("chr", chr_filter)
  protein_ann <- protein_ann[keep_chr, , drop = FALSE]
  if (nrow(protein_ann) == 0) {
    stop("No proteins found for PWAS_CHR_FILTER=", chr_filter, call. = FALSE)
  }
}
if (!is.na(max_proteins) && max_proteins > 0 && nrow(protein_ann) > max_proteins) {
  protein_ann <- protein_ann[seq_len(max_proteins), , drop = FALSE]
}

# Keep numeric covariates only. Exclude IDs and the cell-fraction columns, since
# cell fractions are passed separately through pi_k.
cov_cols <- setdiff(
  names(sample_map),
  c("RawSampleID", "SampleID", "DonorID", "SampleID_pi", "SampleID_cov",
    "Cardiomyocytes", "Other")
)
cov_cols <- cov_cols[vapply(sample_map[cov_cols], is.numeric, logical(1))]
cell_type_cols <- c(
  "Cardiomyocytes",
  "Other"
)

cat("Aligned EA protein samples:", nrow(sample_map), "\n")
cat("Covariates:", paste(cov_cols, collapse = ", "), "\n")
cat("Cell types:", paste(cell_type_cols, collapse = ", "), "\n")
cat("Genotype raw directory:", geno_raw_dir, "\n")
cat("Weight nonzero threshold:", weight_eps, "\n")
cat("MiXcan alpha:", mixcan_alpha, "\n")
cat("Lambda choice:", lambda_choice, "\n")
if (nzchar(chr_filter)) {
  cat("Chromosome filter:", chr_filter, "\n")
}
if (!is.na(max_proteins) && max_proteins > 0) {
  cat("Protein row limit:", max_proteins, "\n")
}

# ------------------------------------------------------------------------------
# Step 5. Train one PWAS MiXcan model per protein
# ------------------------------------------------------------------------------

res_weights_all <- vector("list", nrow(protein_ann))
res_intercepts_all <- vector("list", nrow(protein_ann))
res_tissue_weights_all <- vector("list", nrow(protein_ann))
skipped <- list()
diagnostics <- list()

for (j in seq_len(nrow(protein_ann))) {
  target <- protein_ann[j, ]
  cat("Processing", j, "of", nrow(protein_ann), target$gene_id, "\n")

  # Without gene coordinates we cannot define the cis genotype window.
  if (is.na(target$chr) || is.na(target$start) || is.na(target$end)) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "missing_annotation"
    )
    next
  }

  # Load per-protein or chromosome-wide genotype dosages. Proteins without a
  # matching raw file are logged rather than stopping the full run.
  geno_target <- read_genotype_for_target(
    target_id = target$gene_id,
    chr = target$chr,
    start = target$start,
    end = target$end
  )
  if (is.null(geno_target)) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "missing_or_empty_genotype_raw"
    )
    next
  }

  y_row <- which(protein$gene_id == target$gene_id)
  if (length(y_row) != 1L) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "protein_row_not_unique"
    )
    next
  }

  # Match genotype donors to the aligned protein/pi/covariate donor set.
  donors <- Reduce(
    intersect,
    list(sample_map$DonorID, geno_target$x$DonorID)
  )
  if (length(donors) < 20) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "too_few_genotyped_samples"
    )
    next
  }

  aligned <- sample_map[match(donors, sample_map$DonorID), ]
  geno_rows <- geno_target$x[match(donors, geno_target$x$DonorID), ]

  # Final model matrices:
  #   y    = protein abundance for this protein
  #   x    = cis genotype dosage matrix
  #   z    = covariate matrix
  #   pi_k = two cell-fraction columns: Cardiomyocytes and Other
  y <- as.numeric(protein[y_row, aligned$RawSampleID, drop = TRUE])
  x <- as.matrix(geno_rows[, geno_target$snp_annot$raw_col, drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- donors
  colnames(x) <- geno_target$snp_annot$varID
  candidate_snp_num <- ncol(x)

  z <- as.matrix(aligned[, cov_cols, drop = FALSE])
  pi_k <- as.matrix(aligned[, cell_type_cols, drop = FALSE])

  # Remove samples with any missing value in model inputs.
  complete_idx <- complete.cases(y) &
    complete.cases(x) &
    complete.cases(z) &
    complete.cases(pi_k)
  if (sum(complete_idx) < 20) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "too_few_complete_samples"
    )
    next
  }

  x <- x[complete_idx, , drop = FALSE]
  y <- y[complete_idx]
  z <- z[complete_idx, , drop = FALSE]
  pi_k <- pi_k[complete_idx, , drop = FALSE]
  complete_sample_num <- length(y)

  # Drop extremely rare/unused dosage columns. This mirrors the Breast_BCAC_2CellTypes
  # training script's mean-dosage filter.
  keep_snps <- colMeans(x, na.rm = TRUE) > 0.05
  if (!any(keep_snps)) {
    skipped[[length(skipped) + 1L]] <- data.frame(
      gene_id = target$gene_id,
      reason = "no_snps_after_maf_filter"
    )
    next
  }

  x <- x[, keep_snps, drop = FALSE]
  snp_num_after_maf <- ncol(x)
  snp_annot <- geno_target$snp_annot[match(colnames(x), geno_target$snp_annot$varID), ]
  x_name_matrix <- snp_annot %>%
    transmute(
      varID = varID,
      gene_id = target$gene_id,
      gene_name = target$gene_name,
      chr = chr,
      pos = pos,
      ref_allele = ref_allele,
      eff_allele = eff_allele,
      dosed_allele = dosed_allele
    )

  # Use a deterministic but protein-specific CV split so reruns are reproducible.
  set.seed(1334 + j * 149053)
  foldid <- sample(1:10, length(y), replace = TRUE)

  # Fit the two-component MiXcan model. Errors are logged per protein so one
  # failed fit does not stop training for the remaining proteins.
  fit <- tryCatch(
    MiXcan_train_K_select_lambda(
      y = y,
      x = x,
      pi_k = pi_k,
      cov = z,
      xNameMatrix = x_name_matrix,
      yName = target$gene_id,
      foldid = foldid,
      alpha = mixcan_alpha,
      lambda_choice = lambda_choice
    ),
    error = function(e) {
      skipped[[length(skipped) + 1L]] <<- data.frame(
        gene_id = target$gene_id,
        reason = paste0("MiXcan_train_K_failed: ", conditionMessage(e))
      )
      NULL
    }
  )

  if (is.null(fit)) {
    next
  }

  intercepts <- extract_intercepts(fit, target, cell_type_cols)
  res_intercepts_all[[j]] <- intercepts

  # Save only SNPs with at least one meaningfully nonzero cell-component weight.
  # glmnet can return tiny numerical values around 1e-16, which should not be
  # counted as real selected SNP weights.
  weights <- combine_weights(fit, snp_annot, target, cell_type_cols)
  weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
  weights <- weights[rowSums(abs(weights[, weight_cols, drop = FALSE]) > weight_eps) > 0, ]
  if (nrow(weights)) {
    res_weights_all[[j]] <- weights
  }

  tissue_weights <- combine_tissue_weights(fit, snp_annot, target)
  tissue_weights <- tissue_weights[
    abs(tissue_weights$weight_tissue) > weight_eps,
    ,
    drop = FALSE
  ]
  if (nrow(tissue_weights)) {
    res_tissue_weights_all[[j]] <- tissue_weights
  }

  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    gene_id = target$gene_id,
    gene_name = target$gene_name,
    chr = target$chr,
    sample_num = complete_sample_num,
    candidate_snp_num = candidate_snp_num,
    snp_num_after_maf = snp_num_after_maf,
    selected_snp_num = nrow(weights),
    type = fit$type,
    alpha = mixcan_alpha,
    lambda_choice = lambda_choice,
    lambda_cell = fit$lambda_cell,
    lambda_tissue = fit$lambda_tissue,
    cv_mse_cell = fit$cv_mse_cell,
    cv_mse_tissue = fit$cv_mse_tissue,
    cv_r2_cell = fit$cv_r2_cell,
    cv_r2_tissue = fit$cv_r2_tissue
  )

  if (j %% 100 == 0) {
    cat("Processed", j, "proteins\n")
  }
}

# ------------------------------------------------------------------------------
# Step 6. Save heart protein weights and skipped-protein log
# ------------------------------------------------------------------------------

weights_final <- bind_rows(res_weights_all)
weights_file <- file.path(
  results_dir,
  paste0("weights_heart_protein_cardiomyocytes_other", output_suffix, ".csv")
)
fwrite(weights_final, weights_file)

tissue_weights_final <- bind_rows(res_tissue_weights_all)
tissue_weights_file <- file.path(
  results_dir,
  paste0("predixcan_tissue_weights_heart_protein", output_suffix, ".csv")
)
fwrite(tissue_weights_final, tissue_weights_file)

intercepts_final <- bind_rows(res_intercepts_all)
intercepts_file <- file.path(
  results_dir,
  paste0("intercepts_heart_protein_cardiomyocytes_other", output_suffix, ".csv")
)
fwrite(intercepts_final, intercepts_file)

weight_cols <- paste0("weight_", make_clean_names(cell_type_cols))
snp_terms_final <- weights_final %>%
  mutate(term_type = "snp") %>%
  select(
    gene_id, gene_name, term_type, varID, chr, pos, ref_allele, eff_allele,
    dosed_allele, all_of(weight_cols), type
  )
model_terms_final <- bind_rows(
  intercepts_as_model_terms(intercepts_final, cell_type_cols),
  snp_terms_final
)
model_terms_file <- file.path(
  results_dir,
  paste0("model_terms_heart_protein_cardiomyocytes_other", output_suffix, ".csv")
)
fwrite(model_terms_final, model_terms_file)

skipped_final <- bind_rows(skipped)
skipped_file <- file.path(
  results_dir,
  paste0("weights_heart_protein_cardiomyocytes_other_skipped", output_suffix, ".csv")
)
fwrite(skipped_final, skipped_file)

diagnostics_final <- bind_rows(diagnostics)
diagnostics_file <- file.path(
  results_dir,
  paste0("weights_heart_protein_cardiomyocytes_other_diagnostics", output_suffix, ".csv")
)
fwrite(diagnostics_final, diagnostics_file)

cat("Saved weights:", weights_file, "\n")
cat("Saved PrediXcan tissue weights:", tissue_weights_file, "\n")
cat("Saved intercepts:", intercepts_file, "\n")
cat("Saved model terms:", model_terms_file, "\n")
cat("Saved skipped log:", skipped_file, "\n")
cat("Saved diagnostics:", diagnostics_file, "\n")

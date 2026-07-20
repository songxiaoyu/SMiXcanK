# --- Stage 1: prepare_mw_gwas_input.R ---

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Set directories for gwas data and output
paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
results_dir <- file.path(paper_dir, "Results")
gwas_dir <- file.path(paper_dir, "Data", "BCAC")
output_dir <- file.path(results_dir, "3pi_workspace")

# Read MiXcan model weights
mw_input <- fread(file.path(results_dir, "weights_miXcan_full_pi3.csv"))

# Parse varID to 1000Genome format and create 'CHR' column
mw_input <- mw_input %>%
  mutate(
    CHR = str_extract(varID, "chr\\d+"),
    varID = varID %>%
      str_replace("^chr", "") %>%
      str_replace("_b38$", "") %>%
      str_replace_all("_", ":")
  )
setkey(mw_input, CHR, varID)


# Process each chromosome
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")

  gwas_file <- file.path(gwas_dir, sprintf("chr%d_icogs_hg38.csv", chr))
  #gwas_file <- file.path(gwas_dir, sprintf("chrX_icogs_hg38.csv", chr))
  if (!file.exists(gwas_file)) {
    cat("  --> File not found, skipping\n")
    next
  }
  gwas_df <- fread(gwas_file)[, c(2:7, 9, 48), with = FALSE][, POS := POS_hg38][, POS_hg38 := NULL]

  gwas_fwd <- copy(gwas_df)
  gwas_fwd[, varID := paste0(sub("^chr","", chr.iCOGs), ":", POS, ":", Baseline.Gwas, ":", Effect.Gwas)]
    gwas_fwd[, flip := FALSE]

  gwas_rev <- copy(gwas_df)
  gwas_rev[, varID := paste0(sub("^chr","", chr.iCOGs), ":", POS, ":", Effect.Gwas, ":", Baseline.Gwas)]
  gwas_rev[, flip := TRUE]

  gwas_long <- rbindlist(list(gwas_fwd, gwas_rev), use.names=TRUE, fill=TRUE)
  mw_gwas_input <- mw_input[CHR == paste0("chr", chr)][gwas_long, on = .(varID), nomatch = 0L]

  mw_gwas_input[flip == TRUE, beta.Gwas := -beta.Gwas]
  mw_gwas_input[flip == TRUE, c("Baseline.Gwas","Effect.Gwas") := list(Effect.Gwas, Baseline.Gwas)]


  # Save varid list
  write.table(data.frame(mw_gwas_input$varID),
              file = file.path(output_dir, sprintf("bcac2020_filtered_id/bcac2020_filtered_chr%d_gwas_id_pi3.txt", chr)),
              col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Keep the existing 3pi workspace folder/file naming so this script matches the
  # archived run artifacts that downstream code already uses.
  saveRDS(mw_gwas_input, file = file.path(output_dir, sprintf("baca2020_input/chr%d_mw_gwas_input_bcac2020_oct.rds", chr)))
}

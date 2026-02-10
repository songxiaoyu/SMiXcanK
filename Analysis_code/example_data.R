#example_data.R
##
## Create minimal built-in example data for the SMiXcanK README:
##   Step 1: pi_estimation_K()
##   Step 2: MiXcan_train_K_symmetric()
##   Step 3: SMiXcan_assoc_test_K()
##   Step 4: primo_pipeline_wrap()
##
## Run once from the package root:
##   source("example_data.R")
##
## This will save objects into /data as .rda via usethis::use_data().



# ---- 0) Minimal bulk expression example (genes x samples) ----
# 6 marker genes (3 per cell type), 8 samples
exprB_example <- matrix(
  c(
    # CellType1 marker genes: GENE_A1..GENE_A3
    8.10, 8.00, 7.90, 8.20, 8.10, 7.95, 8.05, 8.15,  # GENE_A1
    7.50, 7.60, 7.40, 7.70, 7.60, 7.55, 7.45, 7.65,  # GENE_A2
    6.90, 7.00, 6.80, 7.10, 7.00, 6.95, 6.85, 7.05,  # GENE_A3
    # CellType2 marker genes: GENE_B1..GENE_B3
    5.20, 5.10, 5.30, 5.20, 5.10, 5.15, 5.25, 5.05,  # GENE_B1
    4.80, 4.70, 4.90, 4.80, 4.70, 4.75, 4.85, 4.65,  # GENE_B2
    4.30, 4.40, 4.20, 4.50, 4.40, 4.35, 4.25, 4.45   # GENE_B3
  ),
  nrow = 6,
  byrow = TRUE
)

rownames(exprB_example) <- c("GENE_A1", "GENE_A2", "GENE_A3",
                             "GENE_B1", "GENE_B2", "GENE_B3")
colnames(exprB_example) <- paste0("Sample", 1:8)

# ---- 1) Marker list (2–3 marker genes per cell type) ----
markers_example <- list(
  CellType1 = c("GENE_A1", "GENE_A2", "GENE_A3"),
  CellType2 = c("GENE_B1", "GENE_B2", "GENE_B3")
)

# ---- 2) Genotype matrix for model training / LD (samples x SNPs) ----
# 8 samples, 5 SNPs
x_example <- matrix(
  c(
    0, 1, 2, 0, 1,
    1, 0, 1, 2, 0,
    2, 1, 0, 1, 2,
    0, 2, 1, 0, 1,
    1, 1, 2, 1, 0,
    2, 0, 1, 2, 1,
    0, 1, 0, 1, 2,
    1, 2, 1, 0, 0
  ),
  nrow = 8,
  byrow = TRUE
)
rownames(x_example) <- colnames(exprB_example)
colnames(x_example) <- paste0("rs", 1:5)

# ---- 3) One-gene expression vector aligned to samples ----
# Length must match nrow(x_example) / ncol(exprB_example)
y_example <- c(2.10, 2.30, 1.90, 2.20, 2.00, 2.15, 2.05, 2.25)
names(y_example) <- rownames(x_example)

# ---- 4) GWAS summary statistics aligned to SNPs ----
# Length must match ncol(x_example)
gwas_example <- list(
  Beta    = c(0.05, -0.02, 0.03, -0.01, 0.02),
  se_Beta = c(0.01,  0.02, 0.015, 0.02, 0.01)
)
# optional names (helps readers)
names(gwas_example$Beta) <- colnames(x_example)
names(gwas_example$se_Beta) <- colnames(x_example)

# ---- 5) Minimal merged table for Step 4 (PRIMO pipeline) ----
# This is an example of what users would get after genome-wide SMiXcan:
# - marginal p-values per cell type: p_1_ct2, p_2_ct2
# - a joint p-value: p_join_ct2
# - a type label: type_ct2 (CellTypeSpecific / NonSpecific)
merged_example <- data.frame(
  gene = paste0("Gene", 1:12),
  type_ct2 = c(
    "CellTypeSpecific","CellTypeSpecific","CellTypeSpecific","CellTypeSpecific","CellTypeSpecific","CellTypeSpecific",
    "NonSpecific","NonSpecific","NonSpecific","NonSpecific","NonSpecific","NonSpecific"
  ),
  p_1_ct2 = c(1e-4, 0.30, 0.02, 0.60, 0.08, 0.90, 0.70, 0.80, 0.90, 0.40, 0.30, 0.20),
  p_2_ct2 = c(0.40, 1e-5, 0.03, 0.70, 0.50, 0.20, 0.90, 0.70, 0.80, 0.60, 0.40, 0.30),
  p_join_ct2 = c(0.60, 0.50, 0.04, 0.90, 0.30, 0.20, 1e-6, 0.20, 0.60, 0.80, 0.40, 0.10),
  stringsAsFactors = FALSE
)

# ---- 6) Save all example objects into /data ----
usethis::use_data(
  exprB_example,
  markers_example,
  x_example,
  y_example,
  gwas_example,
  merged_example,
  overwrite = TRUE
)

message("Saved example datasets: exprB_example, markers_example, x_example, y_example, gwas_example, merged_example")

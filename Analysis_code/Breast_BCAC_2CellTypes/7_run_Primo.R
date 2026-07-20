library(data.table)
library(dplyr)
library(SMiXcan)
paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
results_dir <- file.path(paper_dir, "Results", "2pi_workspace", "bcac2020_result")
data_dir <- file.path(paper_dir, "Data")
#---Input----
combined_path2 <- file.path(results_dir, "bcac2020_result_pi2.csv")
combined2 <- read.csv(combined_path2)

#---Add anotations for supplementary table (can be skipped if only need Primo)----
ensembl_ref = read.csv(file.path(data_dir, "ensembl38.txt"))
setDT(combined2)
combined2[, gene_id_clean := sub("\\..*$", "", gene_id)]  # drop any ENSG version

# Your Ensembl reference table (the one with Gene.stable.ID and Karyotype.band)
setDT(ensembl_ref)

# Build ENSG -> cytoband map (deduplicate just in case)
ensembl_cyto <- unique(
  ensembl_ref[, .(
    ENSG = sub("\\..*$", "", Gene.stable.ID),
    CYTOBAND = Karyotype.band,
    gene_name = Gene.name
  )],
  by = "ENSG"
)

# Join cytoband by ENSG
combined2 <- ensembl_cyto[combined2, on = .(ENSG = gene_id_clean)]
combined2 = combined2[, ENSG := NULL]

# Put CYTOBAND right after gene_id
setcolorder(combined2, c("gene_name", "gene_id", "CYTOBAND",
                         setdiff(names(combined2), c("gene_name","gene_id","CYTOBAND"))))




#----Run Primo analysis-----
combined2 <- as.data.frame(combined2)
primo2 <- infer_celltype_patterns(
  merged = combined2,
  pvals_names = c("p_1", "p_2"),
  p_join_name = "p_join",
  type_col = "type",
  gene_id_col = "gene_name"
)


primo2$unique_by_cell
primo2$shared_specific_genes # 63
primo2$genes_sig_nonspecific #5
primo2$genes_by_pattern #01: 7; 10: 1




out2 <- primo2$out

# Calculate fwer cutoff
combined2$fwer_p_join <- p.adjust(combined2$p_join, method = "bonferroni")
length(which(combined2$fwer_p_join < 0.05))  #32
0.05 / nrow(combined2) # fwer cutoff

# Calculate FDR cutoff
p <- out2$p_join
fdr_cutoff <- max(p[out2$fdr_p_join < 0.1], na.rm = TRUE)
fdr_cutoff #0.000825722

write.csv(out2, file.path(results_dir, "bcac2020_result_pi2_annotated.csv"), row.names = FALSE)


# Write table S1
out2_rename <- out2 %>%
  dplyr::rename(
    model       = type,
    Z_epi       = Z_1,
    p_epi       = p_1,
    Z_stromal   = Z_2,
    p_stromal   = p_2,
    p_joint     = p_join,
    fdr_p_joint = fdr_p_join,
    prob_00     = post_00,
    prob_10     = post_10,
    prob_01     = post_01,
    prob_11     = post_11
  ) %>%
  dplyr::select(
    gene_name, gene_id, chr, CYTOBAND, model, input_snp_num,
    Z_epi, p_epi, Z_stromal, p_stromal,
    p_joint, fdr_p_joint,
    prob_00, prob_10, prob_01, prob_11
  )


head(out2_rename)
write.csv(out2_rename, file.path(results_dir, "tableS2.csv"), row.names = FALSE)

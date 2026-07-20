library(data.table)
library(dplyr)
library(SMiXcan)
paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
results_dir <- file.path(paper_dir, "Results", "SMiXcanK_results")
data_dir <- file.path(paper_dir, "Data")
#---Input----
combined_path3 <- file.path(results_dir, "bcac2020_result_pi3_02.csv")
combined3 = read.csv(combined_path3)

#---Add anotations----
ensembl_ref = read.csv(file.path(data_dir, "ensembl38.txt"))
setDT(combined3)
combined3[, gene_id_clean := sub("\\..*$", "", gene_id)]  # drop any ENSG version

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
combined3 <- ensembl_cyto[combined3, on = .(ENSG = gene_id_clean)]
combined3 = combined3[, ENSG := NULL]

# Put CYTOBAND right after gene_id
setcolorder(combined3, c("gene_name", "gene_id", "CYTOBAND",
                         setdiff(names(combined3), c("gene_name","gene_id","CYTOBAND"))))



# Calculate fwer cutoff
combined3$fwer_p_join <- p.adjust(combined3$p_join, method = "bonferroni")
length(which(combined3$fwer_p_join < 0.05))  #32
0.05 / nrow(combined3) # fwer cutoff


#----Run Primo analysis-----
combined3 <- as.data.frame(combined3)
primo3 <- infer_celltype_patterns(
  merged = combined3,
  pvals_names = c("p_1", "p_2", "p_3"),
  p_join_name = "p_join",
  type_col = "type",
  gene_id_col = "gene_name"
)


primo3$unique_by_cell
primo3$shared_specific_genes # 12
length(primo3$shared_specific_genes)
primo3$genes_sig_nonspecific #8
primo3$genes_by_pattern #001: 9; 010: 4, 011:1, 110:2, 111:9



out3 <- primo3$out
write.csv(out3, file.path(results_dir, "bcac2020_result_pi3_annotated.csv"), row.names = FALSE)

# Calculate FDR cutoff
p <- out3$p_join
fdr_cutoff <- max(p[primo3$out$fdr_p_join  < 0.1], na.rm=TRUE)
fdr_cutoff #0.000424848


# Write table S1
out3_rename <- out3 %>%
  dplyr::rename(
    model       = type,
    Z_adipocyte = Z_1,
    p_adipocyte = p_1,
    Z_fibroblast= Z_2,
    p_fibroblast = p_2,
    Z_epithelial = Z_3,
    p_epithelial = p_3,
    p_joint     = p_join,
    fdr_p_joint = fdr_p_join
  ) %>%
  dplyr::select(
    gene_name, gene_id, chr, CYTOBAND, model, input_snp_num,
    Z_adipocyte, p_adipocyte, Z_fibroblast, p_fibroblast, Z_epithelial, p_epithelial,
    p_joint, fdr_p_joint,
    post_000, post_100, post_010, post_001, post_110, post_101, post_011, post_111
  )

length(which(out3_rename$fdr_p_joint < 0.1)) #33

head(out3_rename)
write.csv(out3_rename, file.path(results_dir, "tableS3.csv"), row.names = FALSE)

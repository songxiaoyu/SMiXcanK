### 2 cell types
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")

merged$fwer_p_join_ct2 <- p.adjust(merged$p_join_ct2, method = "bonferroni")
length(which(merged$fdr_p_join_ct2 < 0.1))
length(which(merged$fdr_p_join_ct2 < 0.1 & merged$type_ct2 == specific_label))
length(which(merged$fdr_p_join_ct2 < 0.1 & merged$type_ct2 == specific_label & merged$Z_1_ct2 * merged$Z_2_ct2 <0))
length(which(merged$fdr_p_join_ct2 < 0.1 & merged$type_ct2 == specific_label & merged$Z_1_ct2 * merged$Z_2_ct2 >0 & abs(merged$Z_1_ct2) > abs(merged$Z_2_ct2)))
length(which(merged$fdr_p_join_ct2 < 0.1 & merged$type_ct2 == specific_label & merged$Z_1_ct2 * merged$Z_2_ct2 >0 & abs(merged$Z_1_ct2) < abs(merged$Z_2_ct2)))

genes_fwer <- merged[which(merged$fwer_p_join_ct2 < 0.05), 'gene_name']

c("PSMB1", "SLC4A7","CASP8","FAR2","EPN2" ,"EP300","CELF4","CRHR1","GTF3A","PDLIM4","USP37","POC1B","RMC1","WDPCP",
"CD200R1","TRMT61A","STXBP4","CTSW", "SNX32","ZNF169","DPP7",
 "RMI1","RPRML","C8orf33", "FES","PLEKHM1","AC008267.5","DPP3",
"AC073912.1" ,"IQCH-AS1" ,  "AC125494.2" ,"AC016355.1")
# Inputs
fdr_cutoff <- 0.1
specific_label   <- "CellTypeSpecific"
nonspecific_label <- "NonSpecific"
spec_idx_ct2 <- which(
  merged$type_ct2 == 'CellTypeSpecific'
)
merged$fdr_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "BH")
merged$fdr_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "BH")
merged$fdr_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "BH")
merged$fdr_p_join_ct2 <- p.adjust(merged$p_join_ct2, method = "BH")
sig_gene_percentage <- length(
  which(
    merged$fdr_p_1_ct2 < fdr_cutoff |
      merged$fdr_p_2_ct2 < fdr_cutoff
  )
) / nrow(merged)


# 2) Run PRIMO on ALL specific rows (assume all specific)

res_ct2 <- primo_map_classify(
  merged,
  c('p_1_ct2', 'p_2_ct2'),
  type_col = "type_ct2",
  alt_props = sig_gene_percentage/2,
  specific_label = "CellTypeSpecific",
  unspecific_label = "NonSpecific",
)

out_ct2 <- res_ct2$out


# 3) Significant filter (OR rule), applied within each subset
sig_spec_idx <- which(
  merged$type_ct2 == 'CellTypeSpecific' & (merged$fdr_p_1_ct2 < fdr_cutoff | merged$fdr_p_2_ct2 < fdr_cutoff)
)

length(sig_spec_idx)

sig_uns_idx <- which(
  merged$type_ct2 == 'NonSpecific' & merged$fdr_p_join_ct2 < fdr_cutoff
)

length(sig_uns_idx)

# 4) MAP class among significant specific rows
#    Using posterior probs from PRIMO result (same row order as merged_ctspec)
post_sig <- out_ct2[sig_spec_idx,c('V1','V2','V3','V4') , drop = FALSE]
#MAP_class <- max.col(post_sig)

post_non00 <- post_sig[, c("V2","V3","V4"), drop = FALSE]

den <- rowSums(post_non00)
post_cond <- post_non00
ok <- den > 0
post_cond[ok, ] <- sweep(post_non00[ok, , drop = FALSE], 1, den[ok], "/")
post_cond[!ok, ] <- NA


MAP_class_cond <- rep(NA_integer_, nrow(post_cond))
MAP_class_cond[ok] <- max.col(post_cond[ok, , drop = FALSE])


MAP_label_cond <- rep(NA_character_, length(MAP_class_cond))
MAP_label_cond[ok] <- c("10","01","11")[MAP_class_cond[ok]]


# 5) Assign Plot C numbers
n_trait1 <- sum(MAP_label_cond == "10", na.rm = TRUE)
n_trait2 <- sum(MAP_label_cond == "01", na.rm = TRUE)
n_shared_specific <- sum(MAP_label_cond == "11", na.rm = TRUE)


# 6) Nonspecific significant count (your "n_shared_nonspecific")
n_shared_nonspecific <- length(sig_uns_idx)

# 7) Total shared
n_shared_total <- n_shared_specific + n_shared_nonspecific

# (Optional) sanity print
cat("n_trait1 =", n_trait1, "\n")
cat("n_trait2 =", n_trait2, "\n")
cat("n_shared_specific =", n_shared_specific, "\n")
cat("n_shared_nonspecific =", n_shared_nonspecific, "\n")
cat("n_shared_total =", n_shared_total, "\n")



# Inputs
fdr_cutoff <- 0.1
specific_label    <- "CellTypeSpecific"
nonspecific_label <- "NonSpecific"
merged$fwer_p_join_ct3 <- p.adjust(merged$p_join_ct3, method = "bonferroni")
length(which(merged$fwer_p_join_ct3 < 0.05))
merged$fdr_p_1_ct3 <- p.adjust(merged$p_1_ct3, method = "BH")
merged$fdr_p_2_ct3 <- p.adjust(merged$p_2_ct3, method = "BH")
merged$fdr_p_3 <- p.adjust(merged$p_3, method = "BH")
merged$fdr_p_join_ct3 <- p.adjust(merged$p_join_ct3, method = "BH")

length(merged[which(merged$fdr_p_join_ct2 < 0.1 & merged$fdr_p_join_ct3 < 0.1), 'gene_name'])
# 0) 指定 3 个 cell types 的列名
p_cols   <- c("p_1_ct3", "p_2_ct3", "p_3")
fdr_cols <- c("fdr_p_1_ct3", "fdr_p_2_ct3", "fdr_p_3")

# 1) 估计“显著基因比例”（OR 规则）
sig_gene_percentage <- sum(
  merged[[fdr_cols[1]]] < fdr_cutoff |
    merged[[fdr_cols[2]]] < fdr_cutoff |
    merged[[fdr_cols[3]]] < fdr_cutoff,
  na.rm = TRUE
) / nrow(merged)

print(sig_gene_percentage)

# 3) Run PRIMO on ALL rows
# 3 traits/cell types -> 2^3 = 8 posterior patterns in out
res_ct3 <- primo_map_classify(
  merged,
  p_cols,
  type_col = "type_ct3",
  alt_props = sig_gene_percentage/3,   # 简单平摊；如果你有更好的先验也可改
  specific_label = specific_label,
  unspecific_label = nonspecific_label
)

out_ct3 <- res_ct3$out
View(out_ct3)

# 4) Significant filter (OR rule), within each subset
sig_spec_idx <- which(
  merged$type_ct3 == specific_label &
      merged$fdr_p_join_ct3  < fdr_cutoff
)
length(which(merged$fdr_p_join_ct3  < fdr_cutoff))
length(sig_spec_idx)

sig_uns_idx <- which(
  merged$type_ct3 == nonspecific_label &
    merged$fdr_p_join_ct3 < fdr_cutoff
)
length(sig_uns_idx)

# 5) Posterior on significant specific rows
# 3 traits => V1..V8 should exist (8 patterns)
post_sig <- out_ct3[sig_spec_idx, paste0("V", 1:8), drop = FALSE]
post_non000 <- post_sig[, paste0("V", 2:8), drop = FALSE]

den <- rowSums(post_non000)
post_cond <- post_non000
ok <- den > 0
post_cond[ok, ] <- sweep(post_non000[ok, , drop = FALSE], 1, den[ok], "/")
post_cond[!ok, ] <- NA

# 条件 MAP：在非000里找最大
MAP_class_cond <- rep(NA_integer_, nrow(post_cond))
MAP_class_cond[ok] <- max.col(post_cond[ok, , drop = FALSE])  # 1..7 对应 V2..V8
Q <- Primo::make_qmat(1:3)
pattern_labels <- apply(Q, 1, paste0, collapse = "")[2:8]
MAP_label_cond <- rep(NA_character_, length(MAP_class_cond))
MAP_label_cond[ok] <- pattern_labels[MAP_class_cond[ok]]

# 6) # ---- 7 Venn regions for 3 CTs ----
venn_counts <- table(MAP_label_cond)

get_n <- function(key) if (key %in% names(venn_counts)) as.integer(venn_counts[[key]]) else 0L


n_ct1_only <- get_n("100")
n_ct2_only <- get_n("010")
n_ct3_only <- get_n("001")
n_ct1_ct2  <- get_n("110")
n_ct1_ct3  <- get_n("101")
n_ct2_ct3  <- get_n("011")
n_shared_specific <- get_n("111")


# 7) nonspecific significant coun
n_shared_nonspecific <- length(sig_uns_idx)

# 8) Total shared
n_shared_total <- n_shared_specific + n_shared_nonspecific


# 9) 打印
# sanity print
cat("CT1 only:", n_ct1_only, "\n")
cat("CT2 only:", n_ct2_only, "\n")
cat("CT3 only:", n_ct3_only, "\n")
cat("CT1&CT2:", n_ct1_ct2, "\n")
cat("CT1&CT3:", n_ct1_ct3, "\n")
cat("CT2&CT3:", n_ct2_ct3, "\n")
cat("Shared (specific, >=2):", n_shared_specific, "\n")
cat("NonSpecific significant:", n_shared_nonspecific, "\n")

# save result to table S2
pattern_labels0 <- apply(Q, 1, paste0, collapse = "")
primo_labels <- paste0("prob_", pattern_labels0)
colnames(out_ct3)[(ncol(out_ct3)-7):ncol(out_ct3)] <- primo_labels

out_ct3_final <- out_ct3[,  c("gene_name", "gene_id", "chr","CYTOBAND", "type_ct3",
                              "input_snp_num_ct3", "Z_1_ct3", "p_1_ct3", "Z_2_ct3", "p_2_ct3","Z_3", "p_3","p_join_ct3", 'fdr_p_join_ct3', primo_labels)]

write.csv(out_ct3_final, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/tableS2_final.csv', row.names = FALSE)

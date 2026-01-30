library('Primo')
library('SMiXcanK')
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")

# 2 cell types
merged$fwer_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "bonferroni")
merged$fwer_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "bonferroni")
merged$fdr_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "BH")
merged$fdr_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "BH")
merged$fdr_p_1_ct3 <- p.adjust(merged$p_1_ct3, method = "BH")
merged$fdr_p_2_ct3 <- p.adjust(merged$p_2_ct3, method = "BH")
merged$fdr_p_3 <- p.adjust(merged$p_3, method = "BH")

# ==============================================================================
# Numbers for Plot C (computed from your provided logic)
# ==============================================================================

# Inputs
fdr_cutoff <- 0.1
specific_label   <- "CellTypeSpecific"
nonspecific_label <- "NonSpecific"

# 1) Split
merged_ctspec <- merged[merged$type_ct2 == specific_label, , drop = FALSE]
merged_unspec <- merged[merged$type_ct2 == nonspecific_label, , drop = FALSE]

# 2) Run PRIMO on ALL specific rows (assume all specific)
res <- primo_map_all(
  merged_ctspec,
  pvals_names = c("p_1_ct2", "p_2_ct2"),
  alt_props   = c(0.05, 0.05)
)

# 3) Significant filter (OR rule), applied within each subset
sig_spec_idx <- which(
  merged_ctspec$fdr_p_1_ct2 < fdr_cutoff | merged_ctspec$fdr_p_2_ct2 < fdr_cutoff
)

sig_uns_idx <- which(
  merged_unspec$fdr_p_1_ct2 < fdr_cutoff | merged_unspec$fdr_p_2_ct2 < fdr_cutoff
)

# 4) MAP class among significant specific rows
#    Using posterior probs from PRIMO result (same row order as merged_ctspec)
post_sig <- res$primo$post_prob[sig_spec_idx, , drop = FALSE]
MAP_class <- max.col(post_sig)

# 5) Assign Plot C numbers (this matches your earlier interpretation)
#    For 2 traits, PRIMO commonly corresponds to patterns:
#    1=null(00), 2=trait1 only(10), 3=trait2 only(01), 4=both(11)
n_trait1 <- sum(MAP_class == 2, na.rm = TRUE)
n_trait2 <- sum(MAP_class == 3, na.rm = TRUE)
n_shared_specific <- sum(MAP_class == 4, na.rm = TRUE)

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


sub1_ct2 <- cbind(merged_ctspec[, c("gene_name", "gene_id","chr", "CYTOBAND", "type_ct2",
                        "input_snp_num_ct2", "Z_1_ct2", "p_1_ct2", "Z_2_ct2", "p_2_ct2", "p_join_ct2"), ],res$post_prob)


sub2_ct2 <- merged_unspec[, c("gene_name", "gene_id", "chr", "CYTOBAND",  "type_ct2",
                              "input_snp_num_ct2", "Z_1_ct2", "p_1_ct2", "Z_2_ct2", "p_2_ct2", "p_join_ct2"), ]
sub2_ct2[c("1", "2", "3", "4")] <- NA


fs1 <- rbind(sub1_ct2, sub2_ct2)
colnames(fs1)[(ncol(fs1)-3):ncol(fs1)] <- c('prob_00', 'prob_10', 'prob_01', 'prob_11')
colnames(fs1)[5:11] <- c('model',"input_snp_num", "Z_epi", "p_epi","Z_stromal","p_stromal","p_joint")
write.csv(fs1, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/tableS1.csv', row.names = FALSE)


# 3 cell types
# ==============================================================================
# Numbers for Plot C (3 cell types, computed from logic)
# ==============================================================================


# ==============================================================================
# 3 cell types
# Numbers for Plot C (computed from your provided logic)
# ==============================================================================


specific_label    <- "CellTypeSpecific"
nonspecific_label <- "NonSpecific"

# ---- column names (EDIT if your ct3 columns have different names) ----

fdr_cutoff <- 0.1
pvals_names_ct3 <- c("p_1_ct3", "p_2_ct3", "p_3")
fdr_names_ct3   <- c("fdr_p_1_ct2", "fdr_p_2_ct2", "fdr_p_3")

# 1) Split
merged_ctspec <- merged[merged$type_ct3 == specific_label, , drop = FALSE]
merged_unspec <- merged[merged$type_ct3 == nonspecific_label, , drop = FALSE]

# 2) Run PRIMO on ALL specific rows
res <- primo_map_all(
  merged_ctspec,
  pvals_names = pvals_names_ct3,
  alt_props   = rep(0.02, length(pvals_names_ct3))
)

# 3) Significant filter (OR across all 3)
sig_spec_idx <- which(rowSums(merged_ctspec[, fdr_names_ct3, drop = FALSE] < fdr_cutoff, na.rm = TRUE) > 0)
sig_uns_idx  <- which(rowSums(merged_unspec[, fdr_names_ct3, drop = FALSE] < fdr_cutoff, na.rm = TRUE) > 0)

# 4) MAP pattern among significant specific rows
post_sig <- res$primo$post_prob[sig_spec_idx, , drop = FALSE]
MAP_pat  <- max.col(post_sig)  # 1..8

# 5) Convert MAP patterns to 3-bit strings (000..111), then count Venn regions
# pattern_bin rows correspond to PRIMO patterns in order:
# 000, 100, 010, 110, 001, 101, 011, 111  (this is expand.grid default with CT1 varying fastest)
pattern_bin <- as.matrix(expand.grid(
  CT1 = c(0, 1),
  CT2 = c(0, 1),
  CT3 = c(0, 1)
))

pat_str <- apply(pattern_bin[MAP_pat, , drop = FALSE], 1, paste0, collapse = "")

# drop null (000)
pat_str <- pat_str[pat_str != "000"]

venn_counts <- table(pat_str)

get_n <- function(key) if (key %in% names(venn_counts)) as.integer(venn_counts[[key]]) else 0L

# ---- 7 Venn regions for 3 CTs ----
n_ct1_only <- get_n("100")
n_ct2_only <- get_n("010")
n_ct3_only <- get_n("001")

n_ct1_ct2  <- get_n("110")
n_ct1_ct3  <- get_n("101")
n_ct2_ct3  <- get_n("011")

n_ct1_ct2_ct3 <- get_n("111")

# optional summaries
n_any_shared_specific <- n_ct1_ct2 + n_ct1_ct3 + n_ct2_ct3 + n_ct1_ct2_ct3

# 6) NonSpecific significant count (mirror your ct2 logic: one number)
n_shared_nonspecific <- length(sig_uns_idx)

# sanity print
cat("CT1 only:", n_ct1_only, "\n")
cat("CT2 only:", n_ct2_only, "\n")
cat("CT3 only:", n_ct3_only, "\n")
cat("CT1&CT2:", n_ct1_ct2, "\n")
cat("CT1&CT3:", n_ct1_ct3, "\n")
cat("CT2&CT3:", n_ct2_ct3, "\n")
cat("CT1&CT2&CT3:", n_ct1_ct2_ct3, "\n")
cat("Shared (specific, >=2):", n_any_shared_specific, "\n")
cat("NonSpecific significant:", n_shared_nonspecific, "\n")


merged_ctspec3 <- merged[which(merged$type_ct3 == 'CellTypeSpecific'), ]
nrow(merged_ctspec3)
merged_unspec3 <- merged[which(merged$type_ct3 == 'NonSpecific'), ]
pvals3 <- merged_ctspec3[,c('p_1_ct3', 'p_2_ct3', 'p_3')]

nrow(pvals)
res3<- Primo::Primo_pval(
  pvals = pvals3,
  alt_props = c(0.02, 0.02, 0.02))


res_select_fdr3 <-  res3$post_prob[which(merged_ctspec3$fdr_p_1_ct3<0.1 | merged_ctspec3$fdr_p_2_ct3<0.1 | merged_ctspec3$fdr_p_3<0.1), ]
nrow(res_select_fdr3)
MAP_class <- max.col(res_select_fdr3)
table(MAP_class) #### 1  6  5  2  2  3 10
merged_unspec <- merged[which(merged$type_ct3 == 'NonSpecific'), ]
nrow(merged_unspec)

m_s <- merged_unspec[which(merged_unspec$fdr_p_1_ct3<0.1 | merged_unspec$fdr_p_2_ct3<0.1 | merged_unspec3$fdr_p_3<0.1), ]
nrow(m_s) #8



sub1_ct3 <- cbind(merged_ctspec3[, c("gene_name", "gene_id", "chr","CYTOBAND","type_ct3",
                                    "input_snp_num_ct3", "Z_1_ct3", "p_1_ct3", "Z_2_ct3", "p_2_ct3","Z_3", "p_3", "p_join_ct3"), ],res3$post_prob)


sub2_ct3 <- merged_unspec3[, c("gene_name", "gene_id", "chr","CYTOBAND", "type_ct3",
                              "input_snp_num_ct3", "Z_1_ct3", "p_1_ct3", "Z_2_ct3", "p_2_ct3","Z_3", "p_3", "p_join_ct3"), ]
sub2_ct3[c("1", "2", "3", "4","5", "6", "7", "8")] <- NA


fs2 <- rbind(sub1_ct3, sub2_ct3)
patterns <- expand.grid(CT1=c(0,1), CT2=c(0,1), CT3=c(0,1))
primo_labels <- paste0("prob_", patterns$CT1, patterns$CT2, patterns$CT3)
colnames(fs2)[(ncol(fs2)-7):ncol(fs2)] <- primo_labels
write.csv(fs2, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/tableS2.csv', row.names = FALSE)





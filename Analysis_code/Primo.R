library('Primo')
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")

# 2 cell types
merged$fwer_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "bonferroni")
merged$fwer_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "bonferroni")
merged$fdr_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "BH")
merged$fdr_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "BH")
merged$fdr_p_1_ct3 <- p.adjust(merged$p_1_ct3, method = "BH")
merged$fdr_p_2_ct3 <- p.adjust(merged$p_2_ct3, method = "BH")
merged$fdr_p_3 <- p.adjust(merged$p_3, method = "BH")

length(which(merged$fwer_p_ct2 < 0.05))
m1 <- merged[which(merged$fwer_p_1_ct2 < 0.05| merged$fwer_p_2_ct2 < 0.05), ]
ncol(m1)
merged_ctspec <- merged[which(merged$type_ct2 == 'CellTypeSpecific'), ]

pvals <- merged_ctspec[,c('p_1_ct2', 'p_2_ct2')]
#pvals <- merged[,c('p_1_ct2', 'p_2_ct2')]
nrow(pvals)
res<- Primo::Primo_pval(
  pvals = pvals,
  alt_props = c(0.05, 0.05))


res_select_fdr <-  res$post_prob[which(merged_ctspec$fdr_p_1_ct2<0.1 | merged_ctspec$fdr_p_2_ct2<0.1), ]

MAP_class <- max.col(res_select_fdr)
table(MAP_class) ####  result 4 8 59
merged_unspec <- merged[which(merged$type_ct2 == 'NonSpecific'), ]
nrow(merged_unspec)

m_s <- merged_unspec[which(merged_unspec$fdr_p_1_ct2<0.1 | merged_unspec$fdr_p_2_ct2<0.1), ]
nrow(m_s) #5

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





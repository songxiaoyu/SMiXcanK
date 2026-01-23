

#---read data---
combined_path <- '/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/bcac2020_result/bcac2020_result_pi3_02.csv'
combined3 = read.csv(combined_path)

combined_path2 <- '/Users/zhusinan/Downloads/S-MiXcan_code_folder/2pi/bcac2020_result/bcac2020_result_pi2.csv'
combined2 <- read.csv(combined_path2)

#---fwer---
merged <- merge(
  combined3, combined2,
  by = c("gene_name", "gene_id", "chr"),
  suffixes = c("_ct3", "_ct2")
)

merged$fwer_p_1_ct2 <- p.adjust(merged$p_1_ct2, method = "bonferroni")
merged$fwer_p_2_ct2 <- p.adjust(merged$p_2_ct2, method = "bonferroni")

sig_fwer_1 <- merged$gene_name[which(merged$fwer_p_1_ct2 < 0.05)]
length(sig_fwer_1)
sig_fwer_2 <- merged$gene_name[which(merged$fwer_p_2_ct2 < 0.05)]
length(sig_fwer_2)
sig_fwer_join <- merged$gene_name[which(merged$fwer_p_ct2 < 0.05)]
length(sig_fwer_join)


total <- unique(c(sig_fwer_1, sig_fwer_2, sig_fwer_join))
length(total)
merged$fwer_p_ct3 <- p.adjust(merged$p_join_ct3, method = "bonferroni")
merged$fwer_p_ct2 <- p.adjust(merged$p_join_ct2, method = "bonferroni")

colnames(merged)[1:2] <- c('gene_id', 'var_name')

#gene3 = combined3$gene_name[which(combined3$p_join < 0.05)]
#length(gene3)
#gene2 = combined$gene_name[which(combined$p_join < 0.05)]
#length(gene2)
#length(intersect(gene3, gene2))

#---add cytoband----
ensembl_ref = read.csv('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data/ensembl38.txt')

setDT(merged)
merged[, gene_id_clean := sub("\\..*$", "", gene_id)]  # drop any ENSG version

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
merged<- ensembl_cyto[merged, on = .(ENSG = gene_id_clean)]

# Clean up helper col
merged = merged[, ENSG := NULL]

# Put CYTOBAND right after gene_id
setcolorder(merged, c("gene_name", "gene_id", "CYTOBAND",
                           setdiff(names(merged), c("gene_name","gene_id","CYTOBAND"))))
merged <- fread("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")


# new model for prob assignment
probs_posthoc <-fit_celltype_em_2ct(merged$Z_1_ct2,merged$Z_2_ct2,max_iter = 20)
result_probs <- summarize_specificity_2ct(probs_posthoc)
table(result_probs$MAP_class)
probs_posthoc$params$pi

fit_try <- fit_celltype_em_2ct(
  z1 = merged$Z_1_ct2,
  z2 = merged$Z_2_ct2,
  max_iter = 300,
  tol = 1e-7,

  # boundary condition
  fix_rho = 0,

  # 关键：限制 11 不要太“肥”
  tau11_cap = 1.6,

  # 关键：别让 10/01 被压到 0（只是稳定 EM，不是“boost power”）
  alpha = c("00"=200, "10"=80, "01"=80, "11"=1),

  # 初始化也很重要：先给 10/01 一些质量
  init = list(
    pi = c("00"=0.90,"10"=0.04,"01"=0.04,"11"=0.02),
    tau1 = 2.5, tau2 = 2.5,
    tau11_1 = 1.3, tau11_2 = 1.3,
    rho = 0
  ),

  verbose = TRUE
)

res_try <- summarize_specificity_2ct(fit_try)
table(res_try$MAP_class)
fit_try$params$pi
fit_try$params[c("tau1","tau2","tau11_1","tau11_2","rho")]


fit2 <- fit_celltype_em_2ct(
  z1 = merged_1$Z_1_ct2,
  z2 = merged_1$Z_2_ct2,
  fix_rho = 0,
  tau11_cap = 3,
  tau_min = 1.2,     #  1.2 / 1.3
  max_iter = 100,
  verbose = TRUE
)

res2 <- summarize_specificity_2ct(fit2)
table(res2$MAP_class)
fit2$params$pi




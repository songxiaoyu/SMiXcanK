library(SMiXcan)
library(data.table)

dat_alt=fread('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_SMiXcan/Results/simulation/final_bar_summary_heterogeneous/replicate_results/heterogeneous_power_eta4_200rep_full_results.csv')
dat_null=fread('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_SMiXcan/Results/simulation/final_bar_summary_heterogeneous/replicate_results/heterogeneous_type1_2000rep_full_results.csv')
table(dat_alt$scenario_id, dat_alt$scenario)


colnames(dat_alt)
colnames(dat_null)
dat_null$scenario_id="S1"


dat=rbind(dat_alt[,c("p_s_join_1", "p_s_join_2","p_s_join", "s_mode", "scenario_id")],
          dat_null[,c("p_s_join_1", "p_s_join_2","p_s_join", "s_mode", "scenario_id")])

dat=fread('/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_SMiXcan/Results/simulation/power_b0_1_b1_1_b2_1_eta4_heter_200rep/power_fixed_setting_full_results.csv')
table(dat$scenario_id)
table(dat$scenario_id, dat$s_mode)

#dat$type=ifelse(dat$s_mode=="invalid_variance_separate", "NS", "CTS")
dat$gene_ID=as.character(seq(1:nrow(dat)))
prob=infer_celltype_patterns(
  merged=as.data.frame(dat),
  pvals_names=c("p_s_join_1", "p_s_join_2"),
  p_join_name="p_s_join",
  fdr_cutoff = 0.1,
  type_col="s_mode",
  gene_id_col="gene_ID",
  specific_label = "joint",
  unspecific_label = "collinear_separate"
)

res=prob$out
# eta1 = 0.2, eta2 = 0.2
idx=which(res$scenario_id=="S2" & res$s_mode=="joint" & res$p_s_join<0.05)
#res[idx,]
tab=table(res[idx,]$MAP_pattern_nonnull)
tab; tab[3]/sum(tab)

# eta1 = 0.2, eta2 = 0
idx=which(res$scenario_id=="S3" & res$s_mode=="joint" & res$p_s_join<0.05)
#res[idx,]
tab=table(res[idx,]$MAP_pattern_nonnull)
tab; tab[2]/sum(tab)

# eta1 = 0, eta2 = 0.2
idx=which(res$scenario_id=="S4" & res$s_mode=="joint" & res$p_s_join<0.05)
#res[idx,]
tab=table(res[idx,]$MAP_pattern_nonnull)
tab; tab[1]/sum(tab)

# eta1 = -0.2, eta2 = 0.2
idx=which(res$scenario_id=="S5" & res$s_mode=="joint" & res$p_s_join<0.05)
#res[idx,]
tab=table(res[idx,]$MAP_pattern_nonnull)
tab; tab[3]/sum(tab)



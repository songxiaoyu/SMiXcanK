#library(devtools)
#library(BayesDeBulk)
library(dplyr)
library(SMiXcan)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
results_dir <- file.path(paper_dir, "Results")

# Three cell types only  ----
# Create Marker Matrix (from Francesca) ----
markers<-list(NULL)

markers[[1]]<-c("FABP4" ,  "BTNL9" ,  "CD300LG", "GPIHBP1", "INHBB"  , "KIF25" ,  "RBP7"  ,  "SEMA3G" , "TCF15"  , "ADIPOQ" , "TIMP4",   "TNMD")
markers[[2]]<-c("CD36", "PDGFRB", "C5AR2", "S100A4", "CD70", "PDPN", "VIM", "ITGA5", "MME", "PDGFRA", "FAP", "ACTA2")
markers[[3]]<-c("CDH1","KRT1","DEFB4","CAV1","MUC1")

names(markers) <- c("adipose","Fibr","Epi")
res_K <- pi_estimation_K(exprB,
                         markers,
                         seed   = 1,
                         n.iter = 10000,
                         burn.in = 1000)

View(res_K$cell.fraction)
View(res_K$cell.expression)
dat=res_K$cell.fraction%>%  as.data.frame() %>% rownames_to_column("SampleID")
# Previous output name used in older runs: BayesDeBulk_pi_3ct.tsv
write.table(
  dat,
  file = file.path(results_dir, "BayesDeBulk_pi_3ct_GTEx.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  append = FALSE,
  quote = FALSE
)

# Step3: Run MiXcan ------------------------------------------------------------------
library(data.table)
# library(xCell)
library(tidyverse)
library(janitor)
library(SMiXcan)
library(readr)
library(dplyr)
library(glmnet)
library(janitor)
library(tibble)
library(doParallel)
library(dplyr)
library(SMiXcanK)
# Set working path
setwd("/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data")
# Step 1: Load data---------------------
# 1. Load and clean GTEx race data. Select only White people.
gtex_race <- read_csv("gtex_v8_race.csv")
gtex_white <- gtex_race %>% filter(RACE == "White") %>% pull(SUBJID)

# 2. Load data and select breast cancer genes
cov = data.frame(fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
breast = fread("Breast_Mammary_Tissue.v8.normalized_expression.bed") # Extract the expression data only to save space
ensembl38 <-read_csv("ensembl38.txt") %>% clean_names()
ensembl38 = unique(ensembl38[, c("gene_stable_id", "gene_name")])
breast$gene_id2 = matrix(unlist(strsplit(breast$gene_id, '[.]')), ncol = 2, byrow = T)[, 1]
breast2 = merge(ensembl38, breast, by.x = "gene_stable_id", by.y = "gene_id2")
dup = unique(breast2[duplicated(breast2$gene_name), "gene_name"])
breast3 = breast2[-which(breast2$gene_name %in% dup), ]

# 3. Select breast cancer expression level for white female
exprB = breast3[, colnames(breast3) %in% cov[which(cov$SEX == 2), "SUBJID"]] # We only need female data
exprB <- exprB %>% dplyr::select(intersect(names(exprB), gtex_white)) # only white
rownames(exprB) = breast3$gene_id
dim(exprB)

# 4. load pi estimation results
pis<- read.csv("/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Data/BayesDeBulk_3CT_GTEx/BayesDeBulk_pi_3ct_GTEx.tsv",
               sep = "\t", header = TRUE)

pis_new <- pis[, c('SampleID','Epi')]
pis_new$Other <- 1- pis_new$Epi
# 5. load GTEx covariate data

cov1=data.frame(fread("phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"))
cov0=cov1[,c("SUBJID", "AGE")]
cov2=fread("Breast_Mammary_Tissue.v8.covariates.txt")
cov3=t(cov2[,-1])
colnames(cov3)=data.frame(cov2)[,1]
cov3=data.frame(colnames(cov2)[-1], cov3);colnames(cov3)[1]="SUBJID"
cov3=cov3[,c(1:21,67:69)] # top 15 PEER factors
cov=data.frame(merge(cov0, cov3, by="SUBJID"))
cov <- cov%>%
  filter(SUBJID %in% gtex_white) %>%
  filter(sex==2) %>%
  select(!sex)


# 6. load genotype
geno1=fread("shapeit_data_for_predictdb_variants-r2") # 178698    847
geno=geno1[,c(1:9,match(cov[,1], colnames(geno1))), with=F]

# 7. elastic net model output
filename <- "en_Breast_Mammary_Tissue.db"
library(DBI)
sqlite.driver <- dbDriver("SQLite")
ElasticNet <- dbConnect(sqlite.driver, dbname = filename)
dbListTables(ElasticNet)
ENextra <- dbReadTable(ElasticNet,"extra")
ENextra$gene_id <- matrix(unlist(strsplit(ENextra$gene, '[.]')), ncol = 2, byrow = T)[, 1]
ENweights <- dbReadTable(ElasticNet,"weights")
ENweights$gene_id <- matrix(unlist(strsplit(ENweights$gene, '[.]')), ncol = 2, byrow = T)[, 1]
# overlapping genes
genID=genID1=intersect(ENextra$gene, breast3$gene_id)
G = length(genID)
G#G  6443 correct


# STEP 2: Model training -----
# run (S)-MiXcan GTEx training - fix seed.


result <- vector("list", G)
res_weights_all <- vector("list", length = G)

# main loop
#G
for (j in 1:G){
  print(j)
  # Process gene j
  yName=genID[j]
  gene = ENextra[which(ENextra$gene == yName), "genename"]
  xName=ENweights[which(ENweights$gene==yName), "varID"]

  xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
  nName=cov$"SUBJID" # women
  n=length(nName)
  # Match the sample ID and select SNPs in this gene for input data
  yData=t(exprB[which(rownames(exprB)==yName), match(nName, colnames(exprB))])
  xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with=F])
  zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
  piData=pis_new[match(nName, pis_new$SampleID),2:3]

  class(xData)<-"numeric"
  # Ignore NaN
  cp.idx=complete.cases(xData) & complete.cases(yData)

  #For each SNP (column), compute the mean genotype across samples in cp.idx, keep SNPs with mean > 0.05, and return those columns for the selected rows.
  #If px == 1 (only one SNP), it computes the mean across cp.idx and keeps it if the mean > 0.05; otherwise drops it.
  px=ncol(xData)
  if (px>1) {
    xvar0=which(apply(xData[cp.idx,], 2, function(f) mean(f)>0.05))
    x.complete=xData[cp.idx,xvar0]
  }
  if (px==1) {
    xvar0=1*(mean(xData[cp.idx,])>0.05)
    x.complete=matrix(xData[,xvar0])

  }
  if (ncol(x.complete) == 0 ||is.null(nrow(x.complete))) {next}

  #Based on cp.idx(sample ID without NaN)
  px2=ncol(x.complete)
  z.complete=zData[cp.idx,]
  xz.complete=as.matrix(cbind(x.complete, z.complete))
  y.complete=yData[cp.idx]
  pi.complete=piData[cp.idx, ]
  length(y.complete)
  pz=ncol(z.complete)
  # Set Seed 10 fold cross validation
  set.seed(1334 + j*149053)
  foldid= sample(1:10, length(y.complete), replace=T)

  # MiXcan method
  # to fix this using training_prediction_model
  #ft.sym=SMiXcan::train_prediction_model(y.train=y.complete, x.train=x.complete, pi.train=pi.complete,cov=z.complete, xNameMatrix=xName.all[xvar0,], foldid=foldid)
  ft.sym <- tryCatch({
    # Attempt to run the function
    ft.sym  <- MiXcan_train_K_symmetric(
      y = y.complete,
      x = x.complete,
      pi_k = pi.complete,
      cov = z.complete,
      xNameMatrix = xName.all[xvar0,],
      foldid = foldid,
      alpha = 0.5
    )
    w1 <- ft.sym$beta.SNP.by.cell$Cell1
    w2 <- ft.sym$beta.SNP.by.cell$Cell2
    w <- cbind(w1, w2$weight)
    colnames(w)[6:7] <- c('weight_cell_1', 'weight_cell_2')
    w$type = ft.sym$type
    if (nrow(w)) res_weights_all[[j]] <- w
  }, error = function(e) {
    # If an error is caught, print a message and return a placeholder result
    cat("MiXcan_train_K_symmetric failed for this gene. Error:", conditionMessage(e), "\n")

    # Return a dummy list structure (or NULL) to prevent the main script from crashing
    return(list(
      type = "ErrorSkipped",
      weight.matrix = NULL,
      beta.all.models = NULL
    ))
  })

  if (j %% 200 == 0) cat("Processed", j, "genes\n")
}

#  bind & save after loop
weights_final <- bind_rows(res_weights_all)
filtered_weights <- weights_final[
  weights_final$weight_cell_1 != 0 | weights_final$weight_cell_2 != 0,
]
write_csv(filtered_weights, "weights_miXcan_full_pi2.csv")
pi3 <- read_csv('weights_miXcan_full_pi3.csv')
print(dim(weights_final))
head(weights_final, 10)


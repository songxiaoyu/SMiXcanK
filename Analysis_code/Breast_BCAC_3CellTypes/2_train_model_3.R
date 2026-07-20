# Train the 3-cell-type MiXcan model on GTEx breast tissue ---------------------------

# Package dependencies
library(data.table)
library(xCell)
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
paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
data_dir <- file.path(paper_dir, "Data")
gtex_dir <- file.path(data_dir, "GTEx")
results_dir <- file.path(paper_dir, "Results")

# ------------------------------------------------------------------------------
# 1. Sample selection and expression input
# ------------------------------------------------------------------------------

# Load GTEx race metadata and keep White participants only.
gtex_race <- read_csv(file.path(gtex_dir, "gtex_v8_race.csv"))
gtex_white <- gtex_race %>% filter(RACE == "White") %>% pull(SUBJID)

# Load breast tissue expression and keep the original object layout.
cov = data.frame(fread(file.path(gtex_dir, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")))
breast = fread(file.path(gtex_dir, "Breast_Mammary_Tissue.v8.normalized_expression.bed")) # Extract the expression data only to save space


# Legacy gene-name harmonization block retained for reference.
#ensembl38 <-read_csv("ensembl38.txt") %>% clean_names()
#ensembl38 = unique(ensembl38[, c("gene_stable_id", "gene_name")])
#breast$gene_id2 = matrix(unlist(strsplit(breast$gene_id, '[.]')), ncol = 2, byrow = T)[, 1]
#breast2 = merge(ensembl38, breast, by.x = "gene_stable_id", by.y = "gene_id2")
#dup = unique(breast2[duplicated(breast2$gene_name), "gene_name"])
#breast3 = breast2[-which(breast2$gene_name %in% dup), ]

# Keep breast-tissue expression for White female samples only.
exprB = as.data.frame(breast)[, colnames(breast) %in% cov[which(cov$SEX == 2), "SUBJID"]] # We only need female data
#exprB = breast3[, colnames(breast3) %in% cov[which(cov$SEX == 2), "SUBJID"]] # We only need female data
exprB <- exprB %>% dplyr::select(intersect(names(exprB), gtex_white)) # only white
rownames(exprB) = breast$gene_id
dim(exprB)

# Load estimated 3-cell-type proportions.
pis <- read.csv(
  file.path(results_dir, "BayesDeBulk_pi_3ct_GTEx.tsv"),
  sep = "\t",
  header = TRUE
)

#pis_old <- read_csv("pi_GTEx.csv") %>% as.data.frame()

# ------------------------------------------------------------------------------
# 2. Covariates and genotype input
# ------------------------------------------------------------------------------

# Load GTEx covariates and keep the samples used above.

cov1=data.frame(fread(file.path(gtex_dir, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")))
cov0=cov1[,c("SUBJID", "AGE")]
cov2=fread(file.path(gtex_dir, "Breast_Mammary_Tissue.v8.covariates.txt"))
cov3=t(cov2[,-1])
colnames(cov3)=data.frame(cov2)[,1]
cov3=data.frame(colnames(cov2)[-1], cov3);colnames(cov3)[1]="SUBJID"
cov3=cov3[,c(1:21,67:69)] # top 15 PEER factors
cov=data.frame(merge(cov0, cov3, by="SUBJID"))
cov <- cov%>%
  filter(SUBJID %in% gtex_white) %>%
  filter(sex==2) %>%
  select(!sex)


# Load genotype data for the selected samples.
geno1=fread(file.path(gtex_dir, "shapeit_data_for_predictdb_variants-r2")) # 178698    847
geno=geno1[,c(1:9,match(cov[,1], colnames(geno1))), with=F]

# Load the reference elastic-net model used to define cis-SNP sets.
filename <- file.path(gtex_dir, "en_Breast_Mammary_Tissue.db")
library(DBI)
sqlite.driver <- dbDriver("SQLite")
ElasticNet <- dbConnect(sqlite.driver, dbname = filename)
dbListTables(ElasticNet)
ENextra <- dbReadTable(ElasticNet,"extra")
ENextra$gene_id <- matrix(unlist(strsplit(ENextra$gene, '[.]')), ncol = 2, byrow = T)[, 1]
ENweights <- dbReadTable(ElasticNet,"weights")
ENweights$gene_id <- matrix(unlist(strsplit(ENweights$gene, '[.]')), ncol = 2, byrow = T)[, 1]
# overlapping genes
genID=genID1=intersect(ENextra$gene, breast$gene_id)

G = length(genID)
G#G  6443 correct


# ------------------------------------------------------------------------------
# 3. Gene-by-gene model training
# ------------------------------------------------------------------------------

# Run the 3-cell-type MiXcan training loop with a fixed CV seed per gene.


result <- vector("list", G)
res_weights_all <- vector("list", length = G)

# Main loop across genes.
for (j in 1:G){
  print(j)
  # Current gene and its SNP annotation in the elastic-net reference model.
  yName=genID[j]
  gene = ENextra[which(ENextra$gene == yName), "genename"]
  xName=ENweights[which(ENweights$gene==yName), "varID"]

  xName.all=ENweights[which(ENweights$gene==yName), c("gene", "rsid", "varID", "ref_allele", "eff_allele")]
  nName=cov$"SUBJID" # women
  n=length(nName)

  # Align expression, genotype, covariates, and cell-type proportions by sample ID.
  yData=t(exprB[which(rownames(exprB)==yName), match(nName, colnames(exprB))])
  xData=t(geno[match(xName, geno$ID), match(nName, colnames(geno)), with=F])
  zData=cov[match(nName, cov$SUBJID),-1]; zData=zData[,-ncol(zData)]
  piData=pis[match(nName, pis$SampleID),2:4]

  class(xData)<-"numeric"

  # Keep samples with complete genotype and expression data.
  cp.idx=complete.cases(xData) & complete.cases(yData)

  # Keep SNPs with mean genotype > 0.05 among complete samples.
  # Handle the single-SNP case separately so matrix dimensions stay valid.
  px=ncol(xData)
  if (px>1) {
    xvar0=which(apply(xData[cp.idx,], 2, function(f) mean(f)>0.05))
    x.complete=xData[cp.idx,xvar0]
  }
  if (px==1) {
    if (mean(xData[cp.idx, 1]) > 0.05) {
      xvar0 <- 1L
      x.complete <- matrix(xData[cp.idx, 1], ncol = 1)
    } else {
      xvar0 <- integer(0)
      x.complete <- matrix(numeric(0), nrow = sum(cp.idx), ncol = 0)
    }

  }
  if (ncol(x.complete) == 0 ||is.null(nrow(x.complete))) {next}

  # Build the final inputs after filtering to complete samples.
  px2=ncol(x.complete)
  z.complete=zData[cp.idx,]
  xz.complete=as.matrix(cbind(x.complete, z.complete))
  y.complete=yData[cp.idx]
  pi.complete=piData[cp.idx, ]
  length(y.complete)
  pz=ncol(z.complete)

  # Create a reproducible 10-fold CV split for this gene.
  set.seed(1334 + j*149053)
  foldid= sample(1:10, length(y.complete), replace=T)

  # Fit the 3-cell-type MiXcan model and store the resulting SNP weights.
  # The older train_prediction_model call is kept below for reference.
  #ft.sym=SMiXcan::train_prediction_model(y.train=y.complete, x.train=x.complete, pi.train=pi.complete,cov=z.complete, xNameMatrix=xName.all[xvar0,], foldid=foldid)
  ft.sym <- tryCatch({
    # Main model fit.
    ft.sym  <- MiXcan_train_K(
      y = y.complete,
      x = x.complete,
      pi_k = pi.complete,
      cov = z.complete,
      xNameMatrix = xName.all[xvar0,],
      foldid = foldid,
      alpha = 0.2
    )
    w1 <- ft.sym$beta.SNP.by.cell$Cell1
    w2 <- ft.sym$beta.SNP.by.cell$Cell2
    w3 <- ft.sym$beta.SNP.by.cell$Cell3
    w <- cbind(w1, w2$weight, w3$weight)
    colnames(w)[6:8] <- c('weight_cell_1', 'weight_cell_2', 'weight_cell_3')
    w$type = ft.sym$type
    if (nrow(w)) res_weights_all[[j]] <- w
  }, error = function(e) {
    # Skip failed genes without stopping the full training loop.
    cat("MiXcan_train_K failed for this gene. Error:", conditionMessage(e), "\n")

    # Return a placeholder object so the outer loop can continue.
    return(list(
      type = "ErrorSkipped",
      weight.matrix = NULL,
      beta.all.models = NULL
    ))
  })
  # ft.sym = MiXcan_train_K_symmetric_ROBUST(y=y.complete, x=x.complete, pi_k=pi.complete,cov=z.complete, xNameMatrix=xName.all[xvar0,], foldid=foldid)
  # MiXcan::MiXcan(y=y.complete, x=x.complete, pi=pi.complete,cov=z.complete, xNameMatrix=xName.all[xvar0,], foldid=foldid)
  # combine weights for both cells on the same SNP rows

  if (j %% 200 == 0) cat("Processed", j, "genes\n")
}

# ------------------------------------------------------------------------------
# 4. Combine and save gene-level weights
# ------------------------------------------------------------------------------

weights_final <- bind_rows(res_weights_all)
filtered_weights <- weights_final[
  weights_final$weight_cell_1 != 0 | weights_final$weight_cell_2 != 0 | weights_final$weight_cell_3 != 0,
]
# Keep the current shared-workflow filename aligned with the archived Dropbox result.
# Previous exploratory runs used names such as weights_miXcan_full_pi3_alpha02.csv.
write_csv(filtered_weights, file.path(results_dir, "weights_miXcan_full_pi3.csv"))

# Quick inspection of the combined weight table.
print(dim(weights_final))
head(weights_final, 10)

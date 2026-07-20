
library(dplyr)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
drive_dir <- file.path(paper_dir, "Data", "DRIVE")
results_dir <- file.path(paper_dir, "Results")

mixcan_assoc_separate <- function(outcome, cell1, cell2, family = "binomial") {
  df <- data.frame(
    outcome = as.numeric(outcome),
    cell1 = as.numeric(cell1),
    cell2 = as.numeric(cell2)
  )
  df <- stats::na.omit(df)

  if (nrow(df) == 0L) {
    return(list(
      cell1_p = NA_real_,
      cell2_p = NA_real_,
      p_combined = NA_real_
    ))
  }

  preds <- scale(as.matrix(df[, c("cell1", "cell2")]), center = TRUE, scale = TRUE)
  fit1 <- stats::glm(df$outcome ~ preds[, 1], family = family)
  fit2 <- stats::glm(df$outcome ~ preds[, 2], family = family)
  cs1 <- summary(fit1)$coefficients
  cs2 <- summary(fit2)$coefficients

  get_p <- function(cs) {
    if (nrow(cs) < 2L) {
      return(NA_real_)
    }
    z <- cs[2, "Estimate"] / cs[2, "Std. Error"]
    2 * stats::pnorm(abs(z), lower.tail = FALSE)
  }

  p1 <- get_p(cs1)
  p2 <- get_p(cs2)
  list(
    cell1_p = p1,
    cell2_p = p2,
    p_combined = safe_ACAT(c(p1, p2))
  )
}

# This is a script-level compatibility helper for the DRIVE comparison only.
# The current package exposes SMiXcan_assoc_test_K(), which returns the newer
# generic output format (p_join_vec and p_join). The older DRIVE analysis code
# below expects separate and joint p-values with two-cell-specific names, so we
# translate the package output here rather than adding a special-purpose helper
# to the package itself.
smixcan_assoc_for_drive <- function(W, gwas_results, x_g, n0, n1, family = "binomial") {
  cov_x <- stats::var(x_g)
  sig_l <- sqrt(diag(cov_x))
  K <- ncol(W)

  z_sep <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    wk <- W[, k]
    sig2_gk <- drop(t(wk) %*% cov_x %*% wk)
    if (is.finite(sig2_gk) && sig2_gk > 0) {
      num <- sum((wk * gwas_results$Beta) * sig_l / gwas_results$se_Beta)
      z_sep[k] <- num / sqrt(sig2_gk)
    }
  }

  p_sep <- 2 * stats::pnorm(abs(z_sep), lower.tail = FALSE)
  joint <- SMiXcan_assoc_test_K(W, gwas_results, x_g, n0 = n0, n1 = n1, family = family)

  list(
    p_1_sep = p_sep[1],
    p_2_sep = p_sep[2],
    p_sep = safe_ACAT(p_sep),
    p_1_join = joint$p_join_vec[1],
    p_2_join = joint$p_join_vec[2],
    p_join = joint$p_join
  )
}

# STEP 1 READ DATA-------------
# read MiXcan weight data

mw_input_new <- read.csv(file.path(results_dir, "weights_miXcan_full_pi2.csv"))
colnames(mw_input_new)[2] ='genename'
mw_input <- mw_input_new %>%
  filter(weight_cell_1 != 0 | weight_cell_2 != 0) %>%
  mutate(CHR = sub("^(chr[^_]+).*", "\\1", varID)) %>%
  mutate(POS = as.integer(sub("^chr[^_]+_([0-9]+).*", "\\1", varID)))
mw_input <- mw_input[,c("gene_id","genename","CHR","POS","weight_cell_1","weight_cell_2","type")]


# read Drive genome data
coln <- scan(text = readLines(file.path(drive_dir, "oncoarray_dosages_chr21.txt"), 1), what = "", quiet = TRUE)[-1]
drive_genome = read.table(file.path(drive_dir, "oncoarray_dosages_chr21.txt"))
colnames(drive_genome)[2:ncol(drive_genome)] <- coln
colnames(drive_genome)[1] <- 'CHR'

# read Drive pheno data
drive_pheno = read.csv(file.path(drive_dir, "oncoarray-drive.pheno.csv"))

#get intersected snp POS
input_pos = mw_input$POS
drive_pos = drive_genome$POS
pos_select = intersect(input_pos, drive_pos) #1102/873


# cleaned input data: mw_drive_input
mw_drive_input <- merge(mw_input, drive_genome, by = "POS")
colnames(mw_drive_input)[16:ncol(mw_drive_input)] = coln[9:length(coln)]

#STEP 2 RUN GWAS------------------
X_genome = t(mw_drive_input[,16:ncol(mw_drive_input)])
D = drive_pheno[drive_pheno$subject_index %in%  rownames(X_genome), "affection_status" ]
 # Match rows of drive_pheno to X_genome by subject_index
D <- drive_pheno[match(rownames(X_genome), drive_pheno$subject_index), "affection_status" ]
 # run gwas
drive_gwas_result = run_gwas(X_genome, D, stats::binomial(), method_binomial = "glm")
drive_gwas_result = data.frame(drive_gwas_result)
 # combine gwas results with inpit
drive_input_total = cbind(drive_gwas_result, mw_drive_input)




#STEP 3 RUN MIXCAN & S-MIXCAN--------------
#split DATA BY GENE
drive_gene <-unique(mw_drive_input$genename) #78 90
drive_split_df <- split(drive_input_total, drive_input_total$genename)

drive_result = data.frame(matrix(ncol = 12, nrow = length(drive_gene)))
colnames(drive_result) <- c('p_m_sep_1','p_m_sep_2','p_m_sep','p_m_join_1','p_m_join_2','p_m_join','p_s_sep_1','p_s_sep_2','p_s_sep','p_s_join_1','p_s_join_2','p_s_join')
filtered_list <- list()

for (g in 1:length(drive_gene)) {
  gene = drive_gene[g]
  gene_df <- drive_split_df[[gene]]  # Access each dataframe
  print(paste("Processing gene:", gene))
  # Extract weight columns
  W1 <- gene_df[["weight_cell_1"]]
  W2 <- gene_df[["weight_cell_2"]]

  # Filter SNPs with non-zero weights
  selected_snp <- gene_df[which(W1 != 0 | W2 != 0), ]

  # Get genome with snp id
  x_selected = t(selected_snp[,18:ncol(selected_snp)])

  # Calculate predicted gene expression value y
  y_predicted <- data.frame(matrix(ncol = 3, nrow = ncol(selected_snp)-17))
  colnames(y_predicted) = c('cell_0', 'cell_1','cell_2')
  # colnames(y_predicted) = c('cell_1','cell_2')
  y_predicted$cell_0 <- rep(1,nrow(x_selected))
  y_predicted$cell_1 <- (as.matrix(x_selected) %*% as.matrix(selected_snp$weight_cell_1))[,1]
  y_predicted$cell_2 <- (as.matrix(x_selected) %*% as.matrix(selected_snp$weight_cell_2))[,1]
  # scale y1, y2
  y_predicted[,2:3] = scale(y_predicted[,2:3])
  y_predicted$subject_index = rownames(x_selected)

  # run MiXcan association test
  asso_data = merge(y_predicted, drive_pheno, by='subject_index')

  drive_mixcan_sep <- mixcan_assoc_separate(
    outcome = asso_data$affection_status,
    cell1 = asso_data$cell_1,
    cell2 = asso_data$cell_2,
    family = "binomial"
  )
  drive_mixcan_join <- MiXcan_assoc_test(
    outcome = asso_data$affection_status,
    cell1 = asso_data$cell_1,
    cell2 = asso_data$cell_2,
    family = "binomial"
  )

  

  drive_result[g, 'p_m_sep_1'] = drive_mixcan_sep$cell1_p
  drive_result[g, 'p_m_sep_2'] = drive_mixcan_sep$cell2_p
  drive_result[g, 'p_m_sep'] = drive_mixcan_sep$p_combined
  drive_result[g, 'p_m_join_1'] = drive_mixcan_join$cell1_p
  drive_result[g, 'p_m_join_2'] = drive_mixcan_join$cell2_p
  drive_result[g, 'p_m_join'] = drive_mixcan_join$p_combined

  #run S-MiXcan
  gwas_results <- list()
  gwas_results$Beta = selected_snp$Beta
  gwas_results$se_Beta =selected_snp$se_Beta
  W <- cbind(W1, W2)
  S_MiXcan_results <- smixcan_assoc_for_drive(
    W = W,
    gwas_results = gwas_results,
    x_g = x_selected,
    n0 = sum(asso_data$affection_status == 0, na.rm = TRUE),
    n1 = sum(asso_data$affection_status == 1, na.rm = TRUE),
    family = "binomial"
  )

  drive_result[g, 'p_s_sep_1'] = S_MiXcan_results$p_1_sep
  drive_result[g, 'p_s_sep_2'] = S_MiXcan_results$p_2_sep
  drive_result[g, 'p_s_sep'] = S_MiXcan_results$p_sep
  drive_result[g, 'p_s_join_1'] = S_MiXcan_results$p_1_join
  drive_result[g, 'p_s_join_2'] = S_MiXcan_results$p_2_join
  drive_result[g, 'p_s_join'] = S_MiXcan_results$p_join
}

rownames(drive_result) <- drive_gene
write.csv(
  drive_result,
  file.path(results_dir, "drive_result_full_lam_new.csv"),
  row.names = TRUE
)
View(drive_result)
plot(drive_result$p_m_sep, drive_result$p_m_join)
abline(0,1)
plot(drive_result$p_s_sep, drive_result$p_s_join)
abline(0,1)
plot(drive_result$p_m_join, drive_result$p_s_join)
abline(0,1)

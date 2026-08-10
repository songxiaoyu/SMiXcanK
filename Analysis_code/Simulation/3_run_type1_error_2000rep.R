# Run one Type I error setting for MiXcan, S-MiXcan, and PrediXcan.
#
# Default setting:
#   b0 = 1, nonzero b1 = 1, nonzero b2 =  1
#   fixed regularization scale = 0.05
#   eta1 = eta2 = 0, so disease has no expression effect.
#
# Outputs:
#   Results/simulation/type1_b0_1_b1_1_b2_1_heter_2000rep/
#     type1_single_setting_full_results.csv
#     type1_single_setting_summary.csv


# Load functions and set up path ----------
rm(list=ls())
suppressPackageStartupMessages({
  library(data.table)
  library(doRNG)
  library(glmnet)
  library(SMiXcan)
  library(doRNG)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  repo_relative <- file.path(getwd(), "Analysis_code", "Simulation")
  if (dir.exists(repo_relative)) {
    return(normalizePath(repo_relative))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "1_data_generation_function.R"))
source(file.path(script_dir, "2_simu_analysis_function.R"))

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = '/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_SMiXcan'
)
x_pool_path <- Sys.getenv(
  "SMIXCAN_SIM_X_POOL",
  unset = file.path(paper_dir, "Data", "Simulation", "X_pool_filtered.csv")
)
X_pool_filtered <- data.table::fread(x_pool_path, data.table = FALSE)

workflow_name <- Sys.getenv("SMIXCAN_TYPE1_WORKFLOW_NAME", unset = "type1_b0_1_b1_1_b2_1_heter_2000rep")
workflow_dir <- file.path(paper_dir, "Results", "simulation", workflow_name)
dir.create(workflow_dir, recursive = TRUE, showWarnings = FALSE)

# set up parameters ----------
B <- as.integer(Sys.getenv("SMIXCAN_TYPE1_BATCH_SIZE", unset = "200"))
ITR <- as.integer(Sys.getenv("SMIXCAN_TYPE1_ITERATIONS", unset = "10"))
workers <- as.integer(Sys.getenv("SMIXCAN_TYPE1_WORKERS", unset = "10"))
n_train <- as.integer(Sys.getenv("SMIXCAN_TYPE1_N_TRAIN", unset = "300"))
n_test <- as.integer(Sys.getenv("SMIXCAN_TYPE1_N_TEST", unset = "3000"))
nfolds <- as.integer(Sys.getenv("SMIXCAN_TYPE1_NFOLDS", unset = "10"))
alpha_level <- as.numeric(Sys.getenv("SMIXCAN_TYPE1_ALPHA_LEVEL", unset = "0.05"))

b0 <- as.numeric(Sys.getenv("SMIXCAN_TYPE1_B0", unset = "1"))
b1 <- as.numeric(Sys.getenv("SMIXCAN_TYPE1_B1", unset = "1"))
b2 <- as.numeric(Sys.getenv("SMIXCAN_TYPE1_B2", unset = "1"))
group <- Sys.getenv("SMIXCAN_TYPE1_GROUP", unset = "heter xy")
reg_scale <- as.numeric(Sys.getenv("SMIXCAN_TYPE1_REG_SCALE", unset = "0.05"))
snp_start <- as.integer(Sys.getenv("SMIXCAN_TYPE1_SNP_START_COL", unset = "2"))
snp_region <- snp_start:(snp_start + 49L)

result_columns <- c(
  "p_s_sep_1", "p_s_sep_2", "p_s_sep",
  "p_m_sep_1", "p_m_sep_2", "p_m_sep",
  "p_s_join_1", "p_s_join_2", "p_s_join",
  "p_m_join_1", "p_m_join_2", "p_m_join",
  "p_predixcan",
  "m_mode", "s_reg_scale", "s_reg_condition",
  "s_mode"
)


# run analysis and save results  ----------
result <- data.frame(matrix(NA, ncol = length(result_columns), nrow = ITR * B))
colnames(result) <- result_columns

for (i in seq_len(ITR)) {
  message("Generating/running batch ", i, "/", ITR)
  data_batch <- DataGen_newhap_binom(
    mc = B,
    n.train = n_train,
    n.test = n_test,
    p = 50,
    b0 = b0,
    nonzero_beta1 = b1,
    nonzero_beta2 = b2,
    gammas = c(0, 0),
    var1 = 0.25,
    var2 = 0.25,
    seed = i*15947 + 146,
    group = group,
    X_pool = X_pool_filtered,
    snp_region = snp_region
  )

  result <- run_simulation(
    data_batch,
    result,
    i,
    family = "binomial",
    regularization = "fixed",
    reg_scale = reg_scale,
    nfolds = nfolds,
    workers = workers,
    smixcan_ld_reference = "test"
  )
}

result$b0 <- b0
result$nonzero_beta1 <- b1
result$nonzero_beta2 <- b2
result$n_train <- n_train
result$n_test <- n_test
result$group <- group
result$reg_scale <- reg_scale
result$scenario <- "Type I error | eta1 = 0, eta2 = 0"

full_file <- file.path(workflow_dir, "heterogeneous_type1_2000rep_full_results.csv")
data.table::fwrite(result, full_file)


# formatting results -----------

p_cols <- c(
  MiXcan_join = "p_m_join",
  SMiXcan_join = "p_s_join",
  MiXcan_sep = "p_m_sep",
  SMiXcan_sep = "p_s_sep",
  PrediXcan = "p_predixcan"
)
summary_table <- data.table::rbindlist(lapply(names(p_cols), function(method) {
  p <- result[[p_cols[[method]]]]
  n_valid <- sum(is.finite(p))
  n_sig <- sum(is.finite(p) & p < alpha_level)
  data.table::data.table(
    method = method,
    n_replicates = length(p),
    n_valid_p = n_valid,
    valid_p_rate = n_valid / length(p),
    n_significant = n_sig,
    type1_error_all_reps = n_sig / length(p),
    type1_error_valid_p = ifelse(n_valid > 0, n_sig / n_valid, NA_real_),
    alpha = alpha_level,
    b0 = b0,
    nonzero_beta1 = b1,
    nonzero_beta2 = b2,
    n_train = n_train,
    n_test = n_test,
    group = group,
    reg_scale = reg_scale
  )
}), fill = TRUE)

format_counts <- function(x) {
  counts <- table(x, useNA = "ifany")
  paste(names(counts), as.integer(counts), collapse = "; ")
}

mode_summary <- c(
  m_mode_counts = format_counts(result$m_mode),
  s_mode_counts = format_counts(result$s_mode)
)

summary_table$m_mode_counts <- mode_summary[["m_mode_counts"]]
summary_table$s_mode_counts <- mode_summary[["s_mode_counts"]]

summary_file <- file.path(workflow_dir, "type1_single_setting_summary.csv")
data.table::fwrite(summary_table, summary_file)

# Compatibility output used by the ABC summary plotting scripts.
summary_file_compat <- file.path(workflow_dir, "type1_regularization_scale_summary.csv")
data.table::fwrite(summary_table, summary_file_compat)

print(summary_table)
message("Wrote full results: ", full_file)
message("Wrote summary: ", summary_file)
message("Wrote compatibility summary: ", summary_file_compat)

# Reproduce the simulation design and add PrediXcan as a comparator.
#
# The default simulated data model follows the binary disease design used in
# the main simulation figure:
#   y1 = b0 + X b1 + e1
#   y2 =      X b2 + e2
#   y  = pi * y1 + (1 - pi) * y2
#   logit P(D = 1) = eta0 + eta1 * y1 + eta2 * y2
#
# Current default replicate rule:
#   Type I error scenario, eta1 = eta2 = 0: 2000 replicates
#   Power scenarios: 200 replicates
#
# Outputs:
#   Results/simulation/simulation_figure_with_predixcan_latest_reps/condition_results/*.csv
#   Results/simulation/simulation_figure_with_predixcan_latest_reps/simulation_summary_with_predixcan.csv

suppressPackageStartupMessages({
  library(data.table)
  library(doRNG)
  library(glmnet)
  library(SMiXcan)
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
source(file.path(script_dir, "2_run_sim.R"))

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
x_pool_path <- Sys.getenv(
  "SMIXCAN_SIM_X_POOL",
  unset = file.path(paper_dir, "Data", "Simulation", "X_pool_filtered.csv")
)
if (!file.exists(x_pool_path)) {
  stop("Missing haplotype pool CSV. Copy X_pool_filtered.csv to Data/Simulation/ or set SMIXCAN_SIM_X_POOL.")
}
X_pool_filtered <- data.table::fread(x_pool_path, data.table = FALSE)

result_root <- Sys.getenv(
  "SMIXCAN_SIM_RESULT_DIR",
  unset = file.path(paper_dir, "Results", "simulation")
)
workflow_name <- Sys.getenv(
  "SMIXCAN_SIM_WORKFLOW_NAME",
  unset = "simulation_figure_with_predixcan_latest_reps"
)
workflow_dir <- file.path(result_root, workflow_name)
condition_dir <- file.path(workflow_dir, "condition_results")
dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("SMIXCAN_SIM_BATCH_SIZE", unset = "100"))
ITR <- as.integer(Sys.getenv("SMIXCAN_SIM_ITERATIONS", unset = "2"))
type1_B <- as.integer(Sys.getenv("SMIXCAN_SIM_TYPE1_BATCH_SIZE", unset = as.character(B)))
type1_ITR <- as.integer(Sys.getenv("SMIXCAN_SIM_TYPE1_ITERATIONS", unset = "20"))
workers <- as.integer(Sys.getenv("SMIXCAN_SIM_WORKERS", unset = "1"))
nfolds <- as.integer(Sys.getenv("SMIXCAN_SIM_NFOLDS", unset = "10"))
n_test <- as.integer(Sys.getenv("SMIXCAN_SIM_N_TEST", unset = "3000"))
regularization <- Sys.getenv("SMIXCAN_SIM_REGULARIZATION", unset = "fixed")
reg_scale <- as.numeric(Sys.getenv("SMIXCAN_SIM_REG_SCALE", unset = "0.05"))
skip_existing <- tolower(Sys.getenv("SMIXCAN_SIM_SKIP_EXISTING", unset = "true")) %in% c("1", "true", "yes")
sim_family <- "binomial"
smixcan_ld_reference <- match.arg(
  Sys.getenv("SMIXCAN_SIM_LD_REFERENCE", unset = "train"),
  choices = c("train", "test")
)

eta_effect <- as.numeric(Sys.getenv("SMIXCAN_SIM_ETA_EFFECT", unset = "0.2"))
alpha_level <- as.numeric(Sys.getenv("SMIXCAN_SIM_ALPHA_LEVEL", unset = "0.05"))
snp_start <- as.integer(Sys.getenv("SMIXCAN_SIM_SNP_START_COL", unset = "2"))
snp_region <- snp_start:(snp_start + 49L)

exp_disease <- data.frame(
  scenario_id = paste0("S", 1:5),
  scenario = c(
    "Type I error | eta1 = 0, eta2 = 0",
    sprintf("Power | eta1 = %s, eta2 = %s", eta_effect, eta_effect),
    sprintf("Power | eta1 = %s, eta2 = 0", eta_effect),
    sprintf("Power | eta1 = 0, eta2 = %s", eta_effect),
    sprintf("Power | eta1 = -%s, eta2 = %s", eta_effect, eta_effect)
  ),
  eta1 = c(0, eta_effect, eta_effect, 0, -eta_effect),
  eta2 = c(0, eta_effect, 0, eta_effect, eta_effect),
  stringsAsFactors = FALSE
)
exp_disease$facet_label <- exp_disease$scenario

make_panel_configs <- function() {
  b0_values <- c(-2, -1, 0, 1, 2)
  effect_values <- c(0, 0.5, 1, 1.5, 2)
  ntrain_values <- c(100, 150, 200, 250, 300)

  rbind(
    data.frame(panel = "a", panel_label = "Homogeneous SNP-Exp", vary = "b0",
               x_value = b0_values, b0 = b0_values, beta1 = 1, beta2 = 1,
               n_train = 300, group = "homo xy", stringsAsFactors = FALSE),
    data.frame(panel = "b", panel_label = "Heterogeneous SNP-Exp", vary = "b0",
               x_value = b0_values, b0 = b0_values, beta1 = 1, beta2 = 1,
               n_train = 300, group = "heter xy", stringsAsFactors = FALSE),
    data.frame(panel = "c", panel_label = "Heterogeneous SNP-Exp", vary = "nonzero b1",
               x_value = effect_values, b0 = 1, beta1 = effect_values, beta2 = 1,
               n_train = 300, group = "heter xy", stringsAsFactors = FALSE),
    data.frame(panel = "d", panel_label = "Heterogeneous SNP-Exp", vary = "nonzero b2",
               x_value = effect_values, b0 = 1, beta1 = 1, beta2 = effect_values,
               n_train = 300, group = "heter xy", stringsAsFactors = FALSE),
    data.frame(panel = "e", panel_label = "Homogeneous SNP-Exp", vary = "ntrain",
               x_value = ntrain_values, b0 = 1, beta1 = 1, beta2 = 1,
               n_train = ntrain_values, group = "homo xy", stringsAsFactors = FALSE),
    data.frame(panel = "f", panel_label = "Heterogeneous SNP-Exp", vary = "ntrain",
               x_value = ntrain_values, b0 = 1, beta1 = 1, beta2 = 1,
               n_train = ntrain_values, group = "heter xy", stringsAsFactors = FALSE)
  )
}

result_columns <- c(
  "p_s_sep_1", "p_s_sep_2", "p_s_sep",
  "p_m_sep_1", "p_m_sep_2", "p_m_sep",
  "p_s_join_1", "p_s_join_2", "p_s_join",
  "p_m_join_1", "p_m_join_2", "p_m_join",
  "p_predixcan",
  "m_mode", "s_reg_scale", "s_reg_condition", "s_mode"
)

method_map <- c(
  "MiXcan" = "p_m_join",
  "S-MiXcan" = "p_s_join",
  "PrediXcan" = "p_predixcan"
)

safe_name <- function(x) {
  x <- gsub("-", "m", as.character(x), fixed = TRUE)
  x <- gsub("\\.", "p", x)
  x
}

parse_env_list <- function(name, default_values) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default_values)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

run_condition <- function(cfg, scen) {
  condition_B <- if (scen$scenario_id == "S1") type1_B else B
  condition_ITR <- if (scen$scenario_id == "S1") type1_ITR else ITR

  condition_file <- file.path(
    condition_dir,
    paste0(
      "panel_", cfg$panel,
      "_", gsub(" ", "_", cfg$vary),
      "_", safe_name(cfg$x_value),
      "_", scen$scenario_id,
      "_", sim_family,
      "_ld_", smixcan_ld_reference,
      "_eta_", safe_name(eta_effect),
      "_", gsub(" ", "_", cfg$group),
      "_rep_", condition_B * condition_ITR,
      ".csv"
    )
  )

  if (skip_existing && file.exists(condition_file)) {
    return(data.table::fread(condition_file, data.table = FALSE))
  }

  result <- data.frame(matrix(NA, ncol = length(result_columns), nrow = condition_ITR * condition_B))
  colnames(result) <- result_columns

  for (i in seq_len(condition_ITR)) {
    Data_batch <- DataGen_newhap_binom(
      mc = condition_B,
      n.train = cfg$n_train,
      n.test = n_test,
      p = 50,
      b0 = cfg$b0,
      nonzero_beta1 = cfg$beta1,
      nonzero_beta2 = cfg$beta2,
      gammas = c(scen$eta1, scen$eta2),
      var1 = 1,
      var2 = 1,
      seed = i + 1,
      group = cfg$group,
      X_pool = X_pool_filtered,
      snp_region = snp_region
    )

    message(
      "Panel ", cfg$panel,
      " | ", cfg$vary, "=", cfg$x_value,
      " | ", scen$scenario_id,
      " | iteration ", i, "/", condition_ITR,
      " | batch size ", condition_B
    )

    result <- run_simulation(
      Data_batch,
      result,
      i,
      family = sim_family,
      regularization = regularization,
      reg_scale = reg_scale,
      nfolds = nfolds,
      workers = workers,
      smixcan_ld_reference = smixcan_ld_reference
    )
  }

  result$panel <- cfg$panel
  result$panel_label <- cfg$panel_label
  result$vary <- cfg$vary
  result$x_value <- cfg$x_value
  result$b0 <- cfg$b0
  result$nonzero_beta1 <- cfg$beta1
  result$nonzero_beta2 <- cfg$beta2
  result$n_train <- cfg$n_train
  result$group <- cfg$group
  result$scenario_id <- scen$scenario_id
  result$scenario <- scen$facet_label
  result$eta1 <- scen$eta1
  result$eta2 <- scen$eta2
  result$sim_family <- sim_family
  result$eta_effect <- eta_effect
  result$smixcan_ld_reference <- smixcan_ld_reference
  result$target_replicates <- condition_B * condition_ITR

  write.csv(result, condition_file, row.names = FALSE)
  result
}

configs <- make_panel_configs()
panel_keep <- parse_env_list("SMIXCAN_SIM_PANELS", unique(configs$panel))
scenario_keep <- parse_env_list("SMIXCAN_SIM_SCENARIOS", exp_disease$scenario_id)
configs <- configs[configs$panel %in% panel_keep, , drop = FALSE]
exp_disease <- exp_disease[exp_disease$scenario_id %in% scenario_keep, , drop = FALSE]

if (nrow(configs) == 0 || nrow(exp_disease) == 0) {
  stop("No simulation conditions selected. Check SMIXCAN_SIM_PANELS and SMIXCAN_SIM_SCENARIOS.")
}
all_results <- list()
idx <- 1L

for (cfg_id in seq_len(nrow(configs))) {
  cfg <- configs[cfg_id, , drop = FALSE]
  for (scen_id in seq_len(nrow(exp_disease))) {
    scen <- exp_disease[scen_id, , drop = FALSE]
    all_results[[idx]] <- run_condition(cfg, scen)
    idx <- idx + 1L
  }
}

full_results <- data.table::rbindlist(all_results, fill = TRUE)
full_file <- file.path(workflow_dir, "simulation_full_results_with_predixcan.csv")
data.table::fwrite(full_results, full_file)

summary_list <- lapply(names(method_map), function(method) {
  p_col <- method_map[[method]]
  dat <- full_results[, c(
    "panel", "panel_label", "vary", "x_value", "scenario_id", "scenario",
    "b0", "nonzero_beta1", "nonzero_beta2", "n_train", "group",
    "sim_family", "eta_effect", "smixcan_ld_reference", p_col
  ), with = FALSE]
  names(dat)[names(dat) == p_col] <- "p_value"
  dat$method <- method
  dat
})

summary_input <- data.table::rbindlist(summary_list, fill = TRUE)
summary_table <- summary_input[
  ,
  .(
    n_replicates = .N,
    n_valid_p = sum(is.finite(p_value)),
    estimate = sum(is.finite(p_value) & p_value < alpha_level) / .N
  ),
  by = .(panel, panel_label, vary, x_value, scenario_id, scenario,
         b0, nonzero_beta1, nonzero_beta2, n_train, group,
         sim_family, eta_effect, smixcan_ld_reference, method)
]

summary_file <- file.path(workflow_dir, "simulation_summary_with_predixcan.csv")
data.table::fwrite(summary_table, summary_file)

message("Wrote full replicate results: ", full_file)
message("Wrote summary table: ", summary_file)

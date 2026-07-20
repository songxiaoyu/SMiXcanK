# Run power simulation for one fixed SNP-expression setting across four
# expression-disease alternatives.
#
# Fixed setting:
#   b0 = 1, nonzero b1 = 0.5, nonzero b2 = 1
#   n_train = 300, n_test = 3000, heter xy
#
# Four power scenarios:
#   eta1 = 0.2,  eta2 = 0.2
#   eta1 = 0.2,  eta2 = 0
#   eta1 = 0,    eta2 = 0.2
#   eta1 = -0.2, eta2 = 0.2
#
# Outputs:
#   Results/simulation/power_b0_1_b1_0p5_b2_1_eta4_reg0p05_200rep/
#     power_fixed_setting_full_results.csv
#     power_fixed_setting_summary.csv
#     power_fixed_setting_bar_sep_predixcan.{pdf,png}

suppressPackageStartupMessages({
  library(data.table)
  library(doRNG)
  library(glmnet)
  library(SMiXcan)
  library(ggplot2)
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
source(file.path(script_dir, "0_hap_generation.R"))
source(file.path(script_dir, "1_data_generation_binary.R"))
source(file.path(script_dir, "2_run_sim.R"))

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
x_pool_path <- Sys.getenv(
  "SMIXCAN_SIM_X_POOL",
  unset = file.path(paper_dir, "Data", "Simulation", "X_pool_filtered.csv")
)
X_pool_filtered <- data.table::fread(x_pool_path, data.table = FALSE)

workflow_name <- Sys.getenv(
  "SMIXCAN_POWER_FIXED_WORKFLOW_NAME",
  unset = "power_b0_1_b1_0p5_b2_1_eta4_reg0p05_200rep"
)
workflow_dir <- file.path(paper_dir, "Results", "simulation", workflow_name)
dir.create(workflow_dir, recursive = TRUE, showWarnings = FALSE)

B <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_BATCH_SIZE", unset = "100"))
ITR <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_ITERATIONS", unset = "2"))
workers <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_WORKERS", unset = "4"))
n_train <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_N_TRAIN", unset = "300"))
n_test <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_N_TEST", unset = "3000"))
nfolds <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_NFOLDS", unset = "10"))
alpha_level <- as.numeric(Sys.getenv("SMIXCAN_POWER_FIXED_ALPHA_LEVEL", unset = "0.05"))
reg_scale <- as.numeric(Sys.getenv("SMIXCAN_POWER_FIXED_REG_SCALE", unset = "0.05"))

b0 <- as.numeric(Sys.getenv("SMIXCAN_POWER_FIXED_B0", unset = "1"))
b1 <- as.numeric(Sys.getenv("SMIXCAN_POWER_FIXED_B1", unset = "0.5"))
b2 <- as.numeric(Sys.getenv("SMIXCAN_POWER_FIXED_B2", unset = "1"))
group <- Sys.getenv("SMIXCAN_POWER_FIXED_GROUP", unset = "heter xy")
snp_start <- as.integer(Sys.getenv("SMIXCAN_POWER_FIXED_SNP_START_COL", unset = "2"))
snp_region <- snp_start:(snp_start + 49L)

scenarios <- data.table::data.table(
  scenario_id = c("S2", "S3", "S4", "S5"),
  scenario = c(
    "eta1 = 0.2, eta2 = 0.2",
    "eta1 = 0.2, eta2 = 0",
    "eta1 = 0, eta2 = 0.2",
    "eta1 = -0.2, eta2 = 0.2"
  ),
  eta1 = c(0.2, 0.2, 0, -0.2),
  eta2 = c(0.2, 0, 0.2, 0.2)
)

result_columns <- c(
  "p_s_sep_1", "p_s_sep_2", "p_s_sep",
  "p_m_sep_1", "p_m_sep_2", "p_m_sep",
  "p_s_join_1", "p_s_join_2", "p_s_join",
  "p_m_join_1", "p_m_join_2", "p_m_join",
  "p_predixcan",
  "m_mode", "s_reg_scale", "s_reg_condition", "s_mode"
)

all_results <- list()

for (scenario_index in seq_len(nrow(scenarios))) {
  scen <- scenarios[scenario_index]
  message("Running ", scen$scenario_id, ": ", scen$scenario)
  result <- data.frame(matrix(NA, ncol = length(result_columns), nrow = ITR * B))
  colnames(result) <- result_columns

  for (i in seq_len(ITR)) {
    message("  Batch ", i, "/", ITR)
    data_batch <- DataGen_newhap_binom(
      mc = B,
      n.train = n_train,
      n.test = n_test,
      p = 50,
      b0 = b0,
      nonzero_beta1 = b1,
      nonzero_beta2 = b2,
      gammas = c(scen$eta1, scen$eta2),
      var1 = 1,
      var2 = 1,
      seed = 1000L * scenario_index + i,
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

  result$scenario_id <- scen$scenario_id
  result$scenario <- scen$scenario
  result$eta1 <- scen$eta1
  result$eta2 <- scen$eta2
  all_results[[scenario_index]] <- as.data.table(result)
}

full_results <- data.table::rbindlist(all_results, fill = TRUE)
full_results$b0 <- b0
full_results$nonzero_beta1 <- b1
full_results$nonzero_beta2 <- b2
full_results$n_train <- n_train
full_results$n_test <- n_test
full_results$group <- group
full_results$reg_scale <- reg_scale

full_file <- file.path(workflow_dir, "power_fixed_setting_full_results.csv")
data.table::fwrite(full_results, full_file)

p_cols <- c(
  "MiXcan" = "p_m_sep",
  "S-MiXcan" = "p_s_sep",
  "PrediXcan" = "p_predixcan"
)

summary_table <- data.table::rbindlist(lapply(names(p_cols), function(method) {
  p_col <- p_cols[[method]]
  full_results[
    ,
    {
      p <- get(p_col)
      n_valid <- sum(is.finite(p))
      n_sig <- sum(is.finite(p) & p < alpha_level)
      .(
        method = method,
        n_replicates = .N,
        n_valid_p = n_valid,
        valid_p_rate = n_valid / .N,
        n_significant = n_sig,
        power_all_reps = n_sig / .N,
        power_valid_p = ifelse(n_valid > 0, n_sig / n_valid, NA_real_)
      )
    },
    by = .(scenario_id, scenario, eta1, eta2)
  ]
}), fill = TRUE)

summary_table[, `:=`(
  b0 = b0,
  nonzero_beta1 = b1,
  nonzero_beta2 = b2,
  n_train = n_train,
  n_test = n_test,
  group = group,
  reg_scale = reg_scale,
  alpha = alpha_level
)]

summary_file <- file.path(workflow_dir, "power_fixed_setting_summary.csv")
data.table::fwrite(summary_table, summary_file)

power_metric <- Sys.getenv("SMIXCAN_POWER_FIXED_PLOT_METRIC", unset = "power_all_reps")
if (!power_metric %in% c("power_all_reps", "power_valid_p")) {
  stop("SMIXCAN_POWER_FIXED_PLOT_METRIC must be 'power_all_reps' or 'power_valid_p'.")
}
plot_data <- copy(summary_table)
plot_data[, power_to_plot := get(power_metric)]
plot_data$scenario <- factor(plot_data$scenario, levels = scenarios$scenario)
plot_data$method <- factor(plot_data$method, levels = c("MiXcan", "S-MiXcan", "PrediXcan"))

bar_plot <- ggplot(plot_data, aes(x = method, y = power_to_plot, fill = method)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", power_to_plot)), vjust = -0.35, size = 3.5) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_fill_manual(values = c(
    "MiXcan" = "#F06A63",
    "S-MiXcan" = "#22B8B8",
    "PrediXcan" = "#4C78A8"
  )) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "grey70"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "Power",
    fill = "Method",
    title = "Power comparison across expression-disease settings",
    subtitle = sprintf(
      "b0 = %s, b1 = %s, b2 = %s; %s",
      b0,
      b1,
      b2,
      ifelse(power_metric == "power_all_reps", "all-replicate denominator", "valid p-value denominator")
    )
  )

plot_pdf <- file.path(workflow_dir, "power_fixed_setting_bar_sep_predixcan.pdf")
plot_png <- file.path(workflow_dir, "power_fixed_setting_bar_sep_predixcan.png")
ggsave(plot_pdf, bar_plot, width = 12, height = 4.5)
ggsave(plot_png, bar_plot, width = 12, height = 4.5, dpi = 300)

print(summary_table)
message("Wrote full results: ", full_file)
message("Wrote summary: ", summary_file)
message("Wrote plot: ", plot_pdf)
message("Wrote plot: ", plot_png)

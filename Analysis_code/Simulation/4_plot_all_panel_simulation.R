# Plot the simulation figure with PrediXcan added.
#
# Input:
#   Results/simulation/simulation_figure_with_predixcan_latest_reps/
#     simulation_summary_with_predixcan.csv, or
#     simulation_full_results_with_predixcan.csv when plotting separate
#     MiXcan/S-MiXcan p-values.
#
# Output:
#   Results/simulation/simulation_figure_with_predixcan_latest_reps/
#     simulation_figure_with_predixcan*.pdf
#     simulation_figure_with_predixcan*.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
result_root <- Sys.getenv(
  "SMIXCAN_SIM_RESULT_DIR",
  unset = file.path(paper_dir, "Results", "simulation")
)
workflow_name <- Sys.getenv(
  "SMIXCAN_SIM_WORKFLOW_NAME",
  unset = "simulation_figure_with_predixcan_latest_reps"
)
workflow_dir <- file.path(result_root, workflow_name)
summary_file <- Sys.getenv(
  "SMIXCAN_SIM_SUMMARY",
  unset = file.path(workflow_dir, "simulation_summary_with_predixcan.csv")
)
full_file <- Sys.getenv(
  "SMIXCAN_SIM_FULL",
  unset = file.path(workflow_dir, "simulation_full_results_with_predixcan.csv")
)
pvalue_mode <- match.arg(
  Sys.getenv("SMIXCAN_SIM_PVALUE_MODE", unset = "join"),
  choices = c("join", "sep")
)
alpha_level <- as.numeric(Sys.getenv("SMIXCAN_SIM_ALPHA_LEVEL", unset = "0.05"))

if (pvalue_mode == "sep") {
  if (!file.exists(full_file)) {
    stop("Missing full replicate table. Run the simulation driver first.")
  }
  full_results <- data.table::fread(full_file)
  sep_method_map <- c(
    "MiXcan sep" = "p_m_sep",
    "S-MiXcan sep" = "p_s_sep",
    "PrediXcan" = "p_predixcan"
  )
  dat <- data.table::rbindlist(lapply(names(sep_method_map), function(method) {
    p_col <- sep_method_map[[method]]
    x <- full_results[, c(
      "panel", "panel_label", "vary", "x_value", "scenario_id", "scenario",
      "b0", "nonzero_beta1", "nonzero_beta2", "n_train", "group", p_col
    ), with = FALSE]
    names(x)[names(x) == p_col] <- "p_value"
    x$method <- method
    x
  }), fill = TRUE)
  dat <- dat[
    ,
    .(
      n_replicates = .N,
      n_valid_p = sum(is.finite(p_value)),
      estimate = sum(is.finite(p_value) & p_value < alpha_level) / .N
    ),
    by = .(panel, panel_label, vary, x_value, scenario_id, scenario,
           b0, nonzero_beta1, nonzero_beta2, n_train, group, method)
  ]
} else {
  if (!file.exists(summary_file)) {
    stop("Missing summary table. Run the simulation driver first.")
  }
  dat <- data.table::fread(summary_file)
}

panel_levels <- c("a", "b", "c", "d", "e", "f")

# Match the simulation figure layout: homogeneous disease first, then null/type I,
# followed by the two cell-type-specific alternatives and opposite effects.
scenario_order <- c("S2", "S1", "S4", "S3", "S5")
scenario_levels <- dat[
  match(scenario_order, scenario_id),
  scenario
]
scenario_levels <- scenario_levels[!is.na(scenario_levels)]
method_levels <- if (pvalue_mode == "sep") {
  c("MiXcan sep", "S-MiXcan sep", "PrediXcan")
} else {
  c("MiXcan", "S-MiXcan", "PrediXcan")
}

dat$panel <- factor(dat$panel, levels = panel_levels)
dat$scenario <- factor(dat$scenario, levels = scenario_levels)
dat$method <- factor(dat$method, levels = method_levels)
dat$row_label <- paste0(dat$panel, "  ", dat$panel_label, " | ", dat$vary)
row_levels <- unique(dat[order(panel)]$row_label)
dat$row_label <- factor(dat$row_label, levels = row_levels)
dat$facet_label <- factor(
  paste(dat$row_label, dat$scenario, sep = "\n"),
  levels = as.vector(t(outer(row_levels, scenario_levels, paste, sep = "\n")))
)

method_colors <- c(
  "MiXcan" = "#F06A63",
  "S-MiXcan" = "#22B8B8",
  "MiXcan sep" = "#F06A63",
  "S-MiXcan sep" = "#22B8B8",
  "PrediXcan" = "#4C78A8"
)

row_count <- length(unique(dat$row_label))
figure_height <- max(3.2, 2.35 * row_count + 1.1)

# Keep the full simulation-scale y-axis when power reaches high values. For small
# smoke tests where all estimates are near 0, zoom in so the test lines are
# visible instead of being compressed against the x-axis.
y_max <- max(dat$estimate, na.rm = TRUE)
y_upper <- if (is.finite(y_max) && y_max < 0.2) {
  min(0.25, max(0.1, ceiling(y_max * 20) / 20))
} else {
  1
}

base_theme <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_line(linewidth = 0.2, color = "grey92"),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
    strip.background = element_rect(fill = "grey70", color = "grey70"),
    strip.text.x = element_text(color = "white", size = 8),
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 9, face = "bold", hjust = 0)
  )

row_plots <- lapply(seq_along(row_levels), function(i) {
  row_name <- row_levels[[i]]
  row_dat <- dat[row_label == row_name]
  ggplot(row_dat, aes(x = x_value, y = estimate, color = method, group = method)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2.0, alpha = 0.85) +
    facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
    coord_cartesian(ylim = c(-0.01, y_upper)) +
    base_theme +
    labs(
      title = row_name,
      x = if (i == length(row_levels)) "Simulation parameter value" else NULL,
      y = if (i == ceiling(length(row_levels) / 2)) "Type I error / power" else NULL,
      color = "method"
    ) +
    theme(
      legend.position = "top",
      axis.title.x = if (i == length(row_levels)) element_text(size = 9) else element_blank()
    )
})

p <- patchwork::wrap_plots(row_plots, ncol = 1, guides = "collect") &
  theme(legend.position = "top")

file_suffix <- if (pvalue_mode == "sep") "_sep" else ""
pdf_file <- file.path(workflow_dir, paste0("simulation_figure_with_predixcan", file_suffix, ".pdf"))
png_file <- file.path(workflow_dir, paste0("simulation_figure_with_predixcan", file_suffix, ".png"))

ggsave(pdf_file, p, width = 12, height = figure_height)
ggsave(png_file, p, width = 12, height = figure_height, dpi = 300)

message("Wrote figure: ", pdf_file)
message("Wrote figure: ", png_file)

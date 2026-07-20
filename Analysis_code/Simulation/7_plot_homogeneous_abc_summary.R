# Build the same A/B/C bar summary for the homogeneous SNP-Exp setting.
#
# This script does not rerun simulations. It uses the fixed homogeneous
# simulation outputs:
#   b0 = 1, b1 = b2 = 1, ntrain = 300, group = "homo xy"
#
# Panels:
#   A. Valid p-value rate for each method.
#   B. Type I error under eta1 = eta2 = 0.
#   C. Power across four expression-disease settings.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
simulation_dir <- file.path(paper_dir, "Results", "simulation")
type1_dir <- file.path(simulation_dir, "type1_reg_scale_0p05_b0_1_b1_1_b2_1_homo_2000rep")
power_dir <- file.path(simulation_dir, "power_b0_1_b1_1_b2_1_homo_eta4_reg0p05_200rep")
out_dir <- file.path(simulation_dir, "final_bar_summary_homogeneous")
out_summary <- file.path(out_dir, "homogeneous_snp_exp_bar_summary_values.csv")

type1_file <- file.path(type1_dir, "type1_regularization_scale_summary.csv")
power_file <- file.path(power_dir, "power_fixed_setting_summary.csv")

method_levels <- c("PrediXcan", "MiXcan", "S-MiXcan")
method_colors <- c(
  "MiXcan" = "#F06A63",
  "S-MiXcan" = "#22B8B8",
  "PrediXcan" = "#4C78A8"
)

scenario_levels <- c(
  "Power | eta1 = 0.2, eta2 = 0.2",
  "Power | eta1 = 0.2, eta2 = 0",
  "Power | eta1 = 0, eta2 = 0.2",
  "Power | eta1 = -0.2, eta2 = 0.2"
)
scenario_labels <- c(
  "eta1 = 0.2, eta2 = 0.2",
  "eta1 = 0.2, eta2 = 0",
  "eta1 = 0, eta2 = 0.2",
  "eta1 = -0.2, eta2 = 0.2"
)

type1_method_map <- c(
  "MiXcan_sep" = "MiXcan",
  "SMiXcan_sep" = "S-MiXcan",
  "PrediXcan" = "PrediXcan"
)

if (file.exists(type1_file) && file.exists(power_file)) {
  type1_summary <- data.table::fread(type1_file)
  power_summary <- data.table::fread(power_file)

  type1_data <- type1_summary[method %in% names(type1_method_map)]
  type1_data[, method_label := factor(type1_method_map[method], levels = method_levels)]
  type1_data[, estimate := type1_error_valid_p]

  power_data <- power_summary[method %in% method_levels & scenario %in% scenario_labels]
  power_data[, method_label := factor(method, levels = method_levels)]
  power_data[, estimate := power_all_reps]
  power_data[, scenario_label := factor(scenario, levels = scenario_labels)]

  homogeneous_data <- data.table::rbindlist(
    list(
      data.table::copy(type1_data)[, result_type := "type1"],
      data.table::copy(power_data)[, result_type := "power"]
    ),
    fill = TRUE
  )
} else if (file.exists(out_summary)) {
  homogeneous_data <- data.table::fread(out_summary)
  type1_data <- homogeneous_data[result_type == "type1"]
  power_data <- homogeneous_data[result_type == "power"]
  type1_data[, method_label := factor(method_label, levels = method_levels)]
  power_data[, method_label := factor(method_label, levels = method_levels)]
  power_data[, scenario_label := factor(scenario, levels = scenario_labels)]
} else {
  stop("Missing homogeneous ABC source summaries and final summary table.")
}

base_theme <- theme_bw(base_size = 12.5) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "grey70"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12.5)
  )

plot_a <- ggplot(type1_data, aes(x = method_label, y = valid_p_rate, fill = method_label)) +
  geom_col(width = 0.65, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * valid_p_rate)), vjust = -0.35, size = 4.2) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  base_theme +
  labs(
    title = "A. Proportion of proteins with prediction model",
    x = NULL,
    y = "Proportion"
  )

plot_b <- ggplot(type1_data, aes(x = method_label, y = estimate, fill = method_label)) +
  geom_col(width = 0.65, color = "grey25", linewidth = 0.25) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey35", linewidth = 0.55) +
  geom_text(aes(y = pmax(estimate, 0.055), label = sprintf("%.3f", estimate)),
            vjust = -0.25, size = 4.2) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.09), expand = expansion(mult = c(0, 0.06))) +
  base_theme +
  labs(
    title = "B. Type I error",
    x = NULL,
    y = "Type I error"
  )

plot_c <- ggplot(power_data, aes(x = method_label, y = estimate, fill = method_label)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", estimate)), vjust = -0.35, size = 3.6) +
  facet_wrap(~ scenario_label, nrow = 1) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.85), expand = expansion(mult = c(0, 0.08))) +
  base_theme +
  labs(
    title = "C. Power",
    x = NULL,
    y = "Power"
  )

legend_data <- data.table(method_label = factor(method_levels, levels = method_levels), y = 1)
legend_row <- ggplot(legend_data, aes(x = method_label, y = y, color = method_label)) +
  geom_point(alpha = 0, size = 3) +
  scale_color_manual(values = method_colors) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  theme_void(base_size = 12.5) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(color = "Method")

draw_combined <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 3,
      ncol = 1,
      heights = grid::unit(c(0.13, 1, 1.15), "null")
    )
  ))
  print(legend_row, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::pushViewport(grid::viewport(
    layout.pos.row = 2,
    layout.pos.col = 1,
    layout = grid::grid.layout(nrow = 1, ncol = 2)
  ))
  print(plot_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::upViewport()
  print(plot_c, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid::upViewport()
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "homogeneous_snp_exp_bar_summary.pdf")
out_png <- file.path(out_dir, "homogeneous_snp_exp_bar_summary.png")

data.table::fwrite(homogeneous_data, out_summary)
grDevices::pdf(out_pdf, width = 10, height = 5.8)
draw_combined()
grDevices::dev.off()

grDevices::png(out_png, width = 10, height = 5.8, units = "in", res = 300)
draw_combined()
grDevices::dev.off()

message("Wrote homogeneous summary table: ", out_summary)
message("Wrote homogeneous figure: ", out_pdf)
message("Wrote homogeneous figure: ", out_png)

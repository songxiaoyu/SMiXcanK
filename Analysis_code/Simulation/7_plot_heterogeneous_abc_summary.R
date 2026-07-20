# Build the A/B/C bar summary for the heterogeneous SNP-Exp setting.
#
# Panels:
#   A. Valid p-value rate under the null, used as a proxy for how often each
#      method returns a usable prediction/association result.
#   B. Type I error from the 2000-replicate null simulation.
#   C. Power for b1 = 0.5, b2 = 1 across four eta settings.
#
# The power panels use the all-replicate denominator so failed/NA models are
# counted as not detected.

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

type1_dir <- file.path(simulation_dir, "type1_reg_scale_0p05_b0_1_b1_0p5_b2_1_2000rep")
power_0p5_1_dir <- file.path(simulation_dir, "power_b0_1_b1_0p5_b2_1_eta4_reg0p05_200rep")
out_dir <- file.path(simulation_dir, "final_bar_summary_heterogeneous")
out_summary <- file.path(out_dir, "heterogeneous_snp_exp_bar_summary_values.csv")

type1_file <- file.path(type1_dir, "type1_regularization_scale_summary.csv")
power_file <- file.path(power_0p5_1_dir, "power_fixed_setting_summary.csv")

method_levels_raw <- c("PrediXcan", "MiXcan_sep", "SMiXcan_sep")
method_labels <- c("PrediXcan", "MiXcan", "S-MiXcan")
method_colors <- c(
  "MiXcan" = "#F06A63",
  "S-MiXcan" = "#22B8B8",
  "PrediXcan" = "#4C78A8"
)

prep_type1 <- function(dat) {
  out <- dat[method %in% method_levels_raw]
  out[, method_label := factor(method, levels = method_levels_raw, labels = method_labels)]
  out[]
}

if (file.exists(type1_file) && file.exists(power_file)) {
  type1_summary <- data.table::fread(type1_file)
  power_0p5_1 <- data.table::fread(power_file)

  type1_plot_data <- prep_type1(type1_summary)
  type1_plot_data[, `:=`(
    result_type = "type1",
    estimate = type1_error_valid_p,
    scenario = "Type I error | eta1 = 0, eta2 = 0",
    scenario_id = "S1",
    nonzero_beta1 = 0.5,
    nonzero_beta2 = 1,
    group = "heter xy",
    reg_scale = 0.05
  )]
} else if (file.exists(out_summary)) {
  heterogeneous_data <- data.table::fread(out_summary)
  type1_plot_data <- heterogeneous_data[result_type == "type1"]
  power_0p5_1 <- heterogeneous_data[result_type == "power"]
  type1_plot_data[, method_label := factor(method_label, levels = method_labels)]
} else {
  stop("Missing heterogeneous ABC source summaries and final summary table.")
}

base_bar_theme <- theme_bw(base_size = 12.5) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "grey70"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12.5)
  )

plot_a <- ggplot(type1_plot_data, aes(x = method_label, y = valid_p_rate, fill = method_label)) +
  geom_col(width = 0.65, color = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * valid_p_rate)), vjust = -0.35, size = 4.2) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  base_bar_theme +
  labs(
    title = "A. Proportion of proteins with prediction model",
    x = NULL,
    y = "Proportion",
    fill = "Method"
  )

plot_b <- ggplot(type1_plot_data, aes(x = method_label, y = type1_error_valid_p, fill = method_label)) +
  geom_col(width = 0.65, color = "grey25", linewidth = 0.25) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey35", linewidth = 0.55) +
  geom_text(aes(y = pmax(type1_error_valid_p, 0.055), label = sprintf("%.3f", type1_error_valid_p)),
            vjust = -0.25, size = 4.2) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 0.09), expand = expansion(mult = c(0, 0.06))) +
  base_bar_theme +
  labs(
    title = "B. Type I error",
    x = NULL,
    y = "Type I error",
    fill = "Method"
  )

prep_power <- function(dat, setting_label) {
  out <- dat[method %in% method_labels]
  out[, method_label := factor(method, levels = method_labels)]
  out[, scenario := factor(
    scenario,
    levels = c(
      "eta1 = 0.2, eta2 = 0.2",
      "eta1 = 0.2, eta2 = 0",
      "eta1 = 0, eta2 = 0.2",
      "eta1 = -0.2, eta2 = 0.2"
    )
  )]
  out[, setting := setting_label]
  out[]
}

power_plot_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "grey70"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12.5)
  )

make_power_plot <- function(dat, title) {
  ggplot(dat, aes(x = method_label, y = power_all_reps, fill = method_label)) +
    geom_col(width = 0.68, color = "grey25", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.3f", power_all_reps)), vjust = -0.35, size = 3.6) +
    facet_wrap(~ scenario, nrow = 1) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 0.55), expand = expansion(mult = c(0, 0.08))) +
    power_plot_theme +
    labs(
      title = title,
      x = NULL,
      y = "Power",
      fill = "Method"
    )
}

plot_c <- make_power_plot(
  prep_power(power_0p5_1, "b1 = 0.5, b2 = 1"),
  "C. Power"
)

power_summary_final <- prep_power(power_0p5_1, "b1 = 0.5, b2 = 1")
power_summary_final[, `:=`(
  result_type = "power",
  estimate = power_all_reps
)]
heterogeneous_data <- data.table::rbindlist(
  list(type1_plot_data, power_summary_final),
  fill = TRUE
)

legend_data <- data.table(method_label = factor(method_labels, levels = method_labels), y = 1)
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
out_pdf <- file.path(out_dir, "heterogeneous_snp_exp_bar_summary.pdf")
out_png <- file.path(out_dir, "heterogeneous_snp_exp_bar_summary.png")

data.table::fwrite(heterogeneous_data, out_summary)
grDevices::pdf(out_pdf, width = 10, height = 5.8)
draw_combined()
grDevices::dev.off()

grDevices::png(out_png, width = 10, height = 5.8, units = "in", res = 300)
draw_combined()
grDevices::dev.off()

message("Wrote heterogeneous summary table: ", out_summary)
message("Wrote heterogeneous figure: ", out_pdf)
message("Wrote heterogeneous figure: ", out_png)

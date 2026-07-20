# Plot the final HERMES HF PWAS result.
#
# This figure uses the full step 5 table generated from the selected
# regularization scale.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ggforce)
  library(bacon)
})

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hf_pwas",
  "hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
results_dir <- file.path(workspace_dir, "hermes_hf_result")
result_table <- file.path(results_dir, "hermes_hf_table_pwas_fixed_0p200.csv")
figure_dir <- file.path(
  paper_dir,
  "Figure",
  "hermes_hf_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

plot_prefix <- "HERMES_HF_scale_0p2"
base_font <- 12

if (!file.exists(result_table)) {
  stop("Missing HERMES HF result table: ", result_table, call. = FALSE)
}

required_cols <- c(
  "gene_id", "gene_name", "model", "p_joint", "fdr_p_joint",
  "MAP_pattern_nonnull"
)
pwas_annotated <- fread(
  result_table,
  colClasses = list(character = "MAP_pattern_nonnull")
)
missing_cols <- setdiff(required_cols, names(pwas_annotated))
if (length(missing_cols) > 0) {
  stop("Missing columns in result table: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

make_qq_plot <- function(tab) {
  qq_tab <- tab[
    is.finite(p_joint) & !is.na(p_joint) & p_joint > 0 & p_joint < 1 &
      !is.na(gene_name)
  ]
  if (nrow(qq_tab) == 0) {
    stop("No valid p_joint values available for QQ plot.", call. = FALSE)
  }

  setorder(qq_tab, p_joint)

  p_bounded <- pmin(pmax(qq_tab$p_joint, .Machine$double.eps), 1 - .Machine$double.eps)
  z_values <- qnorm(1 - p_bounded, lower.tail = FALSE)
  bc <- NULL
  invisible(capture.output(
    bc <- suppressWarnings(bacon(z_values, na.exclude = TRUE))
  ))
  lambda_val <- inflation(bc)

  df_qq <- data.table(
    obs = -log10(qq_tab$p_joint),
    exp = -log10(ppoints(nrow(qq_tab))),
    gene = as.character(qq_tab$gene_name),
    significant = qq_tab$fdr_p_join < 0.05
  )
  label_df <- df_qq[significant == TRUE]

  ggplot(df_qq, aes(x = exp, y = obs)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
    geom_point(
      data = label_df,
      color = "#4E79A7",
      alpha = 0.95,
      size = 2.4
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      color = "black",
      fontface = "bold",
      size = 3.3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    annotate(
      "text",
      x = 0,
      y = max(df_qq$obs, na.rm = TRUE),
      label = sprintf("lambda[GC] == %.3f", lambda_val),
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 5,
      fontface = "bold"
    ) +
    labs(
      x = expression(bold(Expected ~ -log[10](p))),
      y = expression(bold(Observed ~ -log[10](p)))
    ) +
    theme_classic(base_size = base_font) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}

make_split_venn <- function(tab) {
  n_cardiomyocyte <- sum(tab$MAP_pattern_nonnull == "10", na.rm = TRUE)
  n_other <- sum(tab$MAP_pattern_nonnull == "01", na.rm = TRUE)
  n_shared_specific <- sum(
    tab$MAP_pattern_nonnull == "11" & tab$model == "CellTypeSpecific",
    na.rm = TRUE
  )
  n_shared_nonspecific <- sum(
    tab$fdr_p_joint < 0.1 & tab$model == "NonSpecific",
    na.rm = TRUE
  )

  cat("n_cardiomyocyte =", n_cardiomyocyte, "\n")
  cat("n_other =", n_other, "\n")
  cat("n_shared_specific =", n_shared_specific, "\n")
  cat("n_shared_nonspecific =", n_shared_nonspecific, "\n")
  cat("n_shared_total =", n_shared_specific + n_shared_nonspecific, "\n")

  label_top <- "Cell Type\nSpecific"
  label_bottom <- "Non-specific"

  circles <- data.frame(
    x0 = c(-0.8, 0.8),
    y0 = c(0, 0),
    r = c(2, 2),
    type = c("Cardiomyocyte", "Non-cardiomyocyte")
  )

  ggplot() +
    geom_circle(
      data = circles,
      aes(x0 = x0, y0 = y0, r = r, fill = type, color = type),
      alpha = 0.5,
      linewidth = 0.5
    ) +
    geom_segment(
      aes(x = -0.95, xend = 0.95, y = 0, yend = 0),
      color = "black",
      linewidth = 1.2
    ) +
    annotate("text", x = -1.8, y = 0, label = n_cardiomyocyte,
             size = 5, fontface = "bold", color = "black") +
    annotate("text", x = 1.8, y = 0, label = n_other,
             size = 5, fontface = "bold", color = "black") +
    annotate("text", x = 0, y = 0.7,
             label = paste0(label_top, "\n", n_shared_specific),
             size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
    annotate("text", x = 0, y = -0.7,
             label = paste0(label_bottom, "\n", n_shared_nonspecific),
             size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
    annotate("text", x = -0.8, y = 2.2, label = "Cardiomyocyte",
             size = 5, fontface = "bold", color = "black") +
    annotate("text", x = 0.8, y = 2.2, label = "Non-cardiomyocyte",
             size = 5, fontface = "bold", color = "black") +
    coord_fixed() +
    scale_fill_manual(values = c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")) +
    scale_color_manual(values = c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")) +
    theme_void(base_size = base_font) +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}

plot_qq <- make_qq_plot(pwas_annotated)
plot_venn <- make_split_venn(pwas_annotated)
split_figure <- cowplot::plot_grid(
  plot_qq,
  plot_venn,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1.05, 0.95)
)

split_path <- file.path(figure_dir, paste0(plot_prefix, "_QQ_venn.pdf"))
final_split_path <- file.path(workspace_dir, "FINAL_HERMES_HF_PWAS_QQ_COUNT_VENN.pdf")

ggsave(split_path, split_figure, width = 10.2, height = 5.2)
ggsave(final_split_path, split_figure, width = 10.2, height = 5.2)

cat("Input:", result_table, "\n")
cat("QQ + Venn plot:", split_path, "\n")
cat("Final QQ + Count Venn plot:", final_split_path, "\n")

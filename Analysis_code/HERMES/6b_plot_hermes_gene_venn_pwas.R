# Plot HERMES DCM PWAS significant genes in a gene-name Venn diagram.
#
# This script uses the clean step 5 result table and writes the gene-name
# final figure directly into both the figure folder and the workspace root.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ggforce)
  library(bacon)
})

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_pwas",
  "hermes_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
result_table <- file.path(workspace_dir, "hermes_result", "hermes_table_pwas_fixed_0p005.csv")
figure_dir <- file.path(
  paper_dir,
  "Figure",
  "hermes_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(result_table)) {
  stop("Missing HERMES result table: ", result_table, call. = FALSE)
}

tab <- fread(result_table, colClasses = list(character = "MAP_pattern_nonnull"))
required_cols <- c("gene_name", "model", "p_joint", "fdr_p_joint", "MAP_pattern_nonnull")
missing_cols <- setdiff(required_cols, names(tab))
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
    fdr = qq_tab$fdr_p_joint,
    sig_group = fifelse(
      qq_tab$fdr_p_joint < 0.05,
      "FDR < 0.05",
      fifelse(qq_tab$fdr_p_joint < 0.1, "0.05 <= FDR < 0.1", "Not labelled")
    )
  )
  label_df <- df_qq[sig_group != "Not labelled"]

  ggplot(df_qq, aes(x = exp, y = obs)) +
    geom_abline(intercept = 0, slope = 1, color = "red",
                linetype = "dashed", linewidth = 1) +
    geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
    geom_point(data = label_df, aes(color = sig_group), alpha = 0.95, size = 2.4) +
    geom_text_repel(
      data = label_df,
      aes(label = gene, color = sig_group),
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
    scale_color_manual(
      values = c("FDR < 0.05" = "red", "0.05 <= FDR < 0.1" = "black"),
      breaks = c("FDR < 0.05", "0.05 <= FDR < 0.1")
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}

sig <- tab[is.finite(fdr_p_joint) & fdr_p_joint < 0.1]
cardiomyocyte_genes <- sig[MAP_pattern_nonnull == "10", unique(gene_name)]
other_genes <- sig[MAP_pattern_nonnull == "01", unique(gene_name)]
shared_genes <- sig[
  MAP_pattern_nonnull == "11" |
    (model == "NonSpecific" & is.na(MAP_pattern_nonnull)),
  unique(gene_name)
]

make_gene_labels <- function(genes, tab, x, y, size) {
  genes <- sort(unique(genes))
  if (!length(genes)) {
    return(list())
  }
  gene_tab <- unique(tab[gene_name %in% genes, .(gene_name, fdr_p_joint)])
  gene_tab <- gene_tab[order(fdr_p_joint)]
  gene_tab[, sig_group := fifelse(fdr_p_joint < 0.05, "FDR < 0.05", "0.05 <= FDR < 0.1")]
  gene_tab[, y_pos := y + (seq_len(.N) - mean(seq_len(.N))) * 0.28]

  lapply(seq_len(nrow(gene_tab)), function(i) {
    annotate(
      "text",
      x = x,
      y = gene_tab$y_pos[i],
      label = gene_tab$gene_name[i],
      size = size,
      fontface = "bold",
      color = if (gene_tab$sig_group[i] == "FDR < 0.05") "red" else "black"
    )
  })
}

cat("Cardiomyocyte genes:", paste(cardiomyocyte_genes, collapse = ", "), "\n")
cat("Other genes:", paste(other_genes, collapse = ", "), "\n")
cat("Shared genes:", paste(shared_genes, collapse = ", "), "\n")

circles <- data.frame(
  x0 = c(-1.05, 1.05),
  y0 = c(0, 0),
  r = c(2.15, 2.15),
  type = c("Cardiomyocyte", "Non-cardiomyocyte")
)

plot_gene_venn <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = type),
    alpha = 0.25,
    color = "black",
    linewidth = 0.8
  )

for (layer in make_gene_labels(cardiomyocyte_genes, tab, x = -2.0, y = 0, size = 3.7)) {
  plot_gene_venn <- plot_gene_venn + layer
}
for (layer in make_gene_labels(other_genes, tab, x = 2.0, y = 0, size = 3.7)) {
  plot_gene_venn <- plot_gene_venn + layer
}
for (layer in make_gene_labels(shared_genes, tab, x = 0, y = 0, size = 3.2)) {
  plot_gene_venn <- plot_gene_venn + layer
}

plot_gene_venn <- plot_gene_venn +
  annotate("text", x = -1.75, y = 2.35, label = "Cardiomyocyte",
           size = 4.4, fontface = "bold") +
  annotate("text", x = 1.75, y = 2.35, label = "Non-cardiomyocyte",
           size = 4.4, fontface = "bold") +
  coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-2.55, 2.8)) +
  scale_fill_manual(values = c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

plot_qq <- make_qq_plot(tab)
combined_plot <- cowplot::plot_grid(
  plot_qq,
  plot_gene_venn,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1.0, 1.05)
)

figure_path <- file.path(figure_dir, "HERMES_scale_0p005_QQ_venn_genes.pdf")
final_path <- file.path(workspace_dir, "FINAL_HERMES_PWAS_QQ_GENE_VENN.pdf")
ggsave(figure_path, combined_plot, width = 11.2, height = 5.2)
ggsave(final_path, combined_plot, width = 11.2, height = 5.2)

cat("Input:", result_table, "\n")
cat("QQ + Gene Venn plot:", figure_path, "\n")
cat("Final QQ + Gene Venn plot:", final_path, "\n")

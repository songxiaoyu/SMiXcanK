# Plot HERMES HF PWAS significant genes in a Venn diagram.
#
# This is a gene-name version of the Venn panel. It uses the final full
# step 5 table and places significant genes in Cardiomyocyte, Other, or shared
# regions without splitting shared genes into cell-type-specific/non-specific.

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
  "Results", "hermes_hf_pwas",
  "hermes_hf_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
result_table <- file.path(
  workspace_dir,
  "hermes_hf_result",
  "hermes_hf_table_pwas_fixed_0p200.csv"
)
figure_dir <- file.path(
  paper_dir,
  "Figure",
  "hermes_hf_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(result_table)) {
  stop("Missing HERMES HF result table: ", result_table, call. = FALSE)
}

tab <- fread(result_table, colClasses = list(character = "MAP_pattern_nonnull"))
required_cols <- c("gene_name", "model", "p_joint", "fdr_p_joint", "MAP_pattern_nonnull")
missing_cols <- setdiff(required_cols, names(tab))
if (length(missing_cols) > 0) {
  stop("Missing columns in result table: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

sig <- tab[is.finite(fdr_p_joint) & fdr_p_joint < 0.1]

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
    significant = qq_tab$fdr_p_joint < 0.05
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
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}

cardiomyocyte_genes <- sig[MAP_pattern_nonnull == "10", unique(gene_name)]
other_genes <- sig[MAP_pattern_nonnull == "01", unique(gene_name)]
shared_genes <- sig[
  MAP_pattern_nonnull == "11" |
    (model == "NonSpecific" & is.na(MAP_pattern_nonnull)),
  unique(gene_name)
]

collapse_genes <- function(x) {
  if (!length(x)) {
    return("")
  }
  paste(sort(x), collapse = "\n")
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
  ) +
  annotate(
    "text",
    x = -2.0,
    y = 0,
    label = collapse_genes(cardiomyocyte_genes),
    size = 3.7,
    fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 2.0,
    y = 0,
    label = collapse_genes(other_genes),
    size = 3.7,
    fontface = "bold",
    lineheight = 0.9
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = collapse_genes(shared_genes),
    size = 3.1,
    fontface = "bold",
    lineheight = 0.88
  ) +
  annotate("text", x = -1.75, y = 2.35, label = "Cardiomyocyte", size = 4.4, fontface = "bold") +
  annotate("text", x = 1.75, y = 2.35, label = "Non-cardiomyocyte", size = 4.4, fontface = "bold") +
  coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-2.55, 2.8)) +
  scale_fill_manual(values = c("Cardiomyocyte" = "#4E79A7", "Non-cardiomyocyte" = "#F28E2B")) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

gene_combo_path <- file.path(figure_dir, "HERMES_HF_scale_0p2_QQ_venn_genes.pdf")
final_gene_combo_path <- file.path(workspace_dir, "FINAL_HERMES_HF_PWAS_QQ_GENE_VENN.pdf")

plot_qq <- make_qq_plot(tab)
gene_combo <- cowplot::plot_grid(
  plot_qq,
  plot_gene_venn,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  rel_widths = c(1.0, 1.05)
)

ggsave(gene_combo_path, gene_combo, width = 11.2, height = 5.2)
ggsave(final_gene_combo_path, gene_combo, width = 11.2, height = 5.2)

cat("Input:", result_table, "\n")
cat("QQ + Gene Venn plot:", gene_combo_path, "\n")
cat("Final QQ + Gene Venn plot:", final_gene_combo_path, "\n")

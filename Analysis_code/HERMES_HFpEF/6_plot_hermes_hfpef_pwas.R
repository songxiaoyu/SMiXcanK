# Plot HERMES_HFpEF HFpEF PWAS results.
#
# Outputs:
#   Figure/hermes_hfpef_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_HFpEF_PWAS_fixed_0p005_QQ.pdf
#   Figure/hermes_hfpef_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin/HERMES_HFpEF_PWAS_fixed_0p005_QQ_split_venn.pdf

library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(ggforce)

paper_dir <- "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
workspace_dir <- file.path(
  paper_dir,
  "Results", "hermes_hfpef_pwas",
  "hermes_hfpef_workspace_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
results_dir <- file.path(workspace_dir, "hermes_hfpef_result")
figure_dir <- file.path(
  paper_dir,
  "Figure",
  "hermes_hfpef_pwas_moderate_100kb_r2_0.99_alpha0.5_lambdamin"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

base_font <- 12

plot_prefix <- "HERMES_HFpEF_PWAS_fixed_0p005"
result_path <- file.path(results_dir, "hermes_hfpef_result_pwas.csv")
annotated_path <- file.path(results_dir, "hermes_hfpef_result_pwas_fixed_0p005_annotated.csv")

if (!file.exists(result_path)) {
  stop("Missing HERMES_HFpEF HFpEF PWAS result file: ", result_path, call. = FALSE)
}

pwas_result <- read.csv(result_path)
if (file.exists(annotated_path)) {
  pwas_annotated <- read.csv(annotated_path, colClasses = c(MAP_pattern_nonnull = "character"))
} else {
  pwas_annotated <- NULL
  message("Step 5 annotated Primo result not found. Split Venn panel will be skipped.")
}

# QQ plot.
# Prefer the Primo-annotated result when available so labels and split-Venn
# counts use the same post-processed table. Fall back to raw association
# results when step 5 has not been run yet.
qq_source <- if (!is.null(pwas_annotated)) pwas_annotated else pwas_result
pvals <- qq_source$p_join
gene_names <- qq_source$gene_name
mask <- is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals < 1 &
  !is.na(gene_names)

p_clean <- pvals[mask]
g_clean <- as.character(gene_names[mask])
if (length(p_clean) == 0) {
  stop("No valid p_join values available for QQ plot.", call. = FALSE)
}

p_bounded <- pmin(pmax(p_clean, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(1 - p_bounded, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)

df_qq <- data.frame(
  obs = -log10(sort(p_clean)),
  exp = -log10(ppoints(length(p_clean))),
  gene = g_clean[order(p_clean)]
)
top_genes_df <- head(df_qq, 10)

plot_qq <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red",
              linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
  geom_text_repel(
    data = top_genes_df,
    aes(label = gene),
    fontface = "bold",
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0,
    y = max(df_qq$obs, na.rm = TRUE),
    label = bquote(lambda[GC] == .(round(lambda_val, 3))),
    hjust = 0,
    vjust = 1,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    x = expression(bold(Expected ~ -log[10](p))),
    y = expression(bold(Observed ~ -log[10](p))),
    title = "HERMES_HFpEF HFpEF PWAS"
  ) +
  theme_classic(base_size = base_font) +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

qq_path <- file.path(figure_dir, paste0(plot_prefix, "_QQ.pdf"))
ggsave(qq_path, plot_qq, width = 5, height = 5)

# Primo split Venn. Cell-type-specific counts come from PRIMO MAP patterns;
# non-specific shared counts are FDR-significant NonSpecific models.
if (!is.null(pwas_annotated) && "MAP_pattern_nonnull" %in% names(pwas_annotated)) {
  n_cardiomyocytes <- length(which(pwas_annotated$MAP_pattern_nonnull == "10"))
  n_other <- length(which(pwas_annotated$MAP_pattern_nonnull == "01"))
  n_shared_specific <- length(which(
    pwas_annotated$MAP_pattern_nonnull == "11" &
      pwas_annotated$type == "CellTypeSpecific"
  ))
  n_shared_nonspecific <- length(which(
    pwas_annotated$fdr_p_join < 0.1 &
      pwas_annotated$type == "NonSpecific"
  ))

  cat("n_cardiomyocytes =", n_cardiomyocytes, "\n")
  cat("n_other =", n_other, "\n")
  cat("n_shared_specific =", n_shared_specific, "\n")
  cat("n_shared_nonspecific =", n_shared_nonspecific, "\n")
  cat("n_shared_total =", n_shared_specific + n_shared_nonspecific, "\n")

  circles <- data.frame(
    x0 = c(-0.8, 0.8),
    y0 = c(0, 0),
    r = c(2, 2),
    type = c("Cardiomyocytes", "Other")
  )

  plot_venn <- ggplot() +
    geom_circle(
      data = circles,
      aes(x0 = x0, y0 = y0, r = r, fill = type),
      alpha = 0.25,
      color = "black",
      linewidth = 0.8
    ) +
    annotate("text", x = -1.7, y = 0, label = n_cardiomyocytes, size = 6, fontface = "bold") +
    annotate("text", x = 1.7, y = 0, label = n_other, size = 6, fontface = "bold") +
    annotate("text", x = 0, y = 0.35, label = n_shared_specific, size = 6, fontface = "bold") +
    annotate("text", x = 0, y = -0.35, label = n_shared_nonspecific, size = 6, fontface = "bold") +
    annotate("text", x = -1.4, y = 2.25, label = "Cardiomyocytes", size = 4.5, fontface = "bold") +
    annotate("text", x = 1.4, y = 2.25, label = "Other", size = 4.5, fontface = "bold") +
    coord_fixed(xlim = c(-3.2, 3.2), ylim = c(-2.5, 2.7)) +
    scale_fill_manual(values = c("Cardiomyocytes" = "#4E79A7", "Other" = "#F28E2B")) +
    theme_void(base_size = base_font) +
    theme(legend.position = "none")

  split_figure <- cowplot::plot_grid(
    plot_qq,
    plot_venn,
    ncol = 2,
    labels = c("A", "B"),
    label_size = 16,
    label_fontface = "bold"
  )
  ggsave(
    file.path(figure_dir, paste0(plot_prefix, "_QQ_split_venn.pdf")),
    split_figure,
    width = 10,
    height = 5
  )
}

cat("Result input:", result_path, "\n")
cat("Annotated input:", annotated_path, "\n")
cat("QQ plot:", qq_path, "\n")

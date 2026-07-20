# ==============================================================================
# FULL SCRIPT: Figure 1 = 3 scatter panels (All / Cell type 1 / Cell type 2)
#             Figure 2 = QQ plot (B) + split-venn (C)
# ==============================================================================

# ==============================================================================
# 0. Setup: Libraries & Global Settings
# ==============================================================================
library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(tidyr)
library(Primo)
library(ggforce)

paper_dir <- Sys.getenv(
  "PAPER_SMIXCAN_DIR",
  unset = "/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan"
)
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")

base_font <- 12

# ==============================================================================
# 1. Load Data
# ==============================================================================
out2 <- read.csv(
  file.path(results_dir, "2pi_workspace", "bcac2020_result", "bcac2020_result_pi2_annotated.csv"),
  colClasses = c(MAP_pattern_nonnull = "character")
)

drive <- read.csv(
  file.path(results_dir, "drive_result_full_lam_new.csv"),
  row.names = 1
)

# ==============================================================================
# 2. Figure 1: Scatter panels (All / Cell type 1 / Cell type 2)
# ==============================================================================

make_scatter_panel <- function(df, xcol, ycol, subtitle = NULL, base_font = 12) {
  df_use <- df %>%
    filter(.data[[xcol]] > 0, .data[[xcol]] < 1, !is.na(.data[[xcol]]),
           .data[[ycol]] > 0, .data[[ycol]] < 1, !is.na(.data[[ycol]])) %>%
    mutate(
      logx = -log10(.data[[xcol]]),
      logy = -log10(.data[[ycol]])
    )

  cor_val   <- cor(df_use$logx, df_use$logy, method = "pearson")
  cor_label <- paste0("r = ", round(cor_val, 2))
  axis_limit <- max(c(df_use$logx, df_use$logy), na.rm = TRUE) * 1.05

  ggplot(df_use, aes(x = logx, y = logy)) +
    geom_abline(intercept = 0, slope = 1, color = "red",
                linetype = "dashed", linewidth = 1) +
    geom_point(color = "#4E79A7", alpha = 0.3, size = 1.2) +
    annotate("text", x = 0, y = axis_limit, label = cor_label,
             hjust = 0, vjust = 1, size = 5, fontface = "bold") +
    xlim(0, axis_limit) + ylim(0, axis_limit) +
    labs(
      x = expression(bold(-log[10](p) ~ "w. individual genotype")),
      y = expression(bold(-log[10](p) ~ "w. summary statistics")),
      title = subtitle
    ) +
    theme_classic(base_size = base_font) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_font + 1, hjust = 0.5),
      axis.line  = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Panel A: overall
plot_f1_a <- make_scatter_panel(
  df = drive,
  xcol = "p_m_join",
  ycol = "p_s_join",
  subtitle = "Tissue-level",
  base_font = base_font
)

# Panel B: cell type 1
plot_f1_b <- make_scatter_panel(
  df = drive,
  xcol = "p_m_join_1",
  ycol = "p_s_join_1",
  subtitle = "Epithelial",
  base_font = base_font
)


# Panel C: cell type 2
plot_f1_c <- make_scatter_panel(
  df = drive,
  xcol = "p_m_join_2",
  ycol = "p_s_join_2",
  subtitle = "Stromal",
  base_font = base_font
)

# Align axes for the three panels
f1_aligned <- cowplot::align_plots(
  plot_f1_a, plot_f1_b, plot_f1_c,
  align = "hv",
  axis = "tblr"
)

figure1 <- cowplot::plot_grid(
  f1_aligned[[1]], f1_aligned[[2]], f1_aligned[[3]],
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("A", "B", "C"),
  label_size = 16,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

print(figure1)

ggsave(file.path(figure_dir, "Figure1_ABC_scatter.pdf"),
       figure1, width = 15, height = 5)

# ==============================================================================
# 3. Figure 2 Panel B: QQ Plot (requires combined2 in your environment)
# ==============================================================================

# Use the annotated Primo output directly for the QQ plot.
pvals_b    <- out2$p_join
genename_b <- out2$gene_name
mask_b <- is.finite(pvals_b) & !is.na(pvals_b) & !is.na(genename_b) & pvals_b > 0 & pvals_b < 1
p_clean_b <- pvals_b[mask_b]
g_clean_b <- as.character(genename_b[mask_b])

p <- pmin(pmax(p_clean_b, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(1 - p, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)

df_qq <- data.frame(
  obs  = -log10(sort(p_clean_b)),
  exp  = -log10(ppoints(length(p_clean_b))),
  gene = g_clean_b[order(p_clean_b)]
)
top_genes_df <- head(df_qq, 10)

plot_f2_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
  geom_text_repel(data = top_genes_df, aes(label = gene),
                  fontface = "bold", size = 3.5) +
  annotate("text", x = 0, y = max(df_qq$obs, na.rm = TRUE),
           label = bquote(lambda[GC] == .(round(lambda_val, 3))),
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  labs(
    x = expression(bold(Expected ~ -log[10](p))),
    y = expression(bold(Observed ~ -log[10](p)))
  ) +
  theme_classic(base_size = base_font) +
  theme(
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==============================================================================
# 4. Figure 2 Panel C: Primo numbers + custom split Venn
# ==============================================================================

n_trait1 <- length(which(out2$MAP_pattern_nonnull == "10"))
n_trait2 <- length(which(out2$MAP_pattern_nonnull == "01"))
n_shared_specific <- length(which(out2$MAP_pattern_nonnull == "11" & out2$type == "CellTypeSpecific"))
n_shared_nonspecific <- length(which(out2$fdr_p_join < 0.1 & out2$type == "NonSpecific"))
n_shared_total <- n_shared_specific + n_shared_nonspecific

cat("n_trait1 =", n_trait1, "\n")
cat("n_trait2 =", n_trait2, "\n")
cat("n_shared_specific =", n_shared_specific, "\n")
cat("n_shared_nonspecific =", n_shared_nonspecific, "\n")
cat("n_shared_total =", n_shared_total, "\n")

label_top    <- "Cell Type\nSpecific"
label_bottom <- "Non-specific"

circles <- data.frame(
  x0 = c(-0.8, 0.8),
  y0 = c(0, 0),
  r  = c(2, 2),
  type = c("Epithelial", "Stromal")
)

venn_colors <- c("Epithelial" = "#FDB462", "Stromal" = "#4E79A7")

plot_f2_c <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = type, color = type),
    alpha = 0.5, linewidth = 0.5
  ) +
  scale_fill_manual(values = venn_colors) +
  scale_color_manual(values = venn_colors) +
  geom_segment(aes(x = -0.95, xend = 0.95, y = 0, yend = 0),
               color = "black", linewidth = 1.2) +
  annotate("text", x = -1.8, y = 0, label = paste0(n_trait1),
           size = 5, fontface = "bold", color = "black") +
  annotate("text", x =  1.8, y = 0, label = paste0(n_trait2),
           size = 5, fontface = "bold", color = "black") +
  annotate("text", x = 0, y = 0.7,
           label = paste0(label_top, "\n", n_shared_specific),
           size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
  annotate("text", x = 0, y = -0.7,
           label = paste0(label_bottom, "\n", n_shared_nonspecific),
           size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
  annotate("text", x = -0.8, y = 2.2, label = "Epithelial",
           size = 5, fontface = "bold", color = "black") +
  annotate("text", x =  0.8, y = 2.2, label = "Stromal",
           size = 5, fontface = "bold", color = "black") +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# Make C a bit taller to match the QQ plot height
plot_f2_c_tall <- cowplot::plot_grid(
  ggdraw(), plot_f2_c, ggdraw(),
  ncol = 1,
  rel_heights = c(0.18, 1, 0.18)
)

# ==============================================================================
# 5. Assemble Figure 2: B + C
# ==============================================================================

bc_aligned <- cowplot::align_plots(
  plot_f2_b, plot_f2_c_tall,
  align = "hv",
  axis  = "tblr"
)

figure2 <- cowplot::plot_grid(
  bc_aligned[[1]], bc_aligned[[2]],
  ncol = 2,
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

print(figure2)

ggsave(file.path(figure_dir, "Figure2_BC.pdf"),
       figure2, width = 10, height = 5)


# Assemble Figure - ABC combination ver
# ==============================================================================
# 5. Final assembly: merge Figure 1 and Figure 2 into one figure
#    Top row:    A (contains the 3 scatter subplots)
#    Bottom row: B (QQ plot) and C (split-Venn)
# ==============================================================================

final_figure <- cowplot::plot_grid(
  figure1,
  cowplot::plot_grid(
    plot_f2_b, plot_f2_c_tall,
    ncol = 2,
    rel_widths = c(1, 1),
    labels = c("B", "C"),
    label_size = 16,
    label_fontface = "bold",
    label_x = 0, label_y = 1
  ),
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("A", ""),
  label_size = 16,
  label_fontface = "bold",
  label_x = 0, label_y = 1
)

print(final_figure)

ggsave(
  file.path(figure_dir, "Final_Figure_ABC.pdf"),
  final_figure,
  width = 15,
  height = 10
)

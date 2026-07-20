# ==============================================================================
# UPDATE: Figure 1 becomes 4 panels (All + ct1 + ct2 + ct3), like your first script
#   Figure 1 = 4 scatter panels (All / Cell type 1 / Cell type 2 / Cell type 3)
#   Figure 2 = QQ (A) + Venn (B)   [unchanged from prior reply]
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(cowplot)
library(bacon)
library(dplyr)
library(tidyr)
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
out3 <- read.csv(
  file.path(results_dir, "SMiXcanK_results", "bcac2020_result_pi3_annotated.csv"),
  colClasses = c(MAP_pattern_nonnull = "character")
)
drive <- read.csv(
  file.path(results_dir, "drive_result_pi3.csv")
)

# ==============================================================================
# 2. Helper: scatter panel
# ==============================================================================
make_scatter_panel <- function(df, xcol, ycol, subtitle = NULL, base_font = 12) {
  if (!all(c(xcol, ycol) %in% names(df))) {
    stop(sprintf("Missing columns in drive: %s / %s", xcol, ycol))
  }

  df_use <- df %>%
    filter(.data[[xcol]] > 0, .data[[xcol]] < 1, !is.na(.data[[xcol]]),
           .data[[ycol]] > 0, .data[[ycol]] < 1, !is.na(.data[[ycol]])) %>%
    mutate(
      logx = -log10(.data[[xcol]]),
      logy = -log10(.data[[ycol]])
    )

  if (nrow(df_use) < 3) {
    cor_label  <- "r = NA"
    axis_limit <- 1
  } else {
    cor_val   <- cor(df_use$logx, df_use$logy, method = "pearson")
    cor_label <- paste0("r = ", round(cor_val, 2))
    axis_limit <- max(c(df_use$logx, df_use$logy), na.rm = TRUE) * 1.05
  }

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

# ==============================================================================
# 3. Figure 1: 4 scatter panels (All / ct1 / ct2 / ct3)
# ==============================================================================
plot_f1_all <- make_scatter_panel(drive, "p_m_join",   "p_s_join",   subtitle = "Tissue-level",         base_font = base_font)
plot_f1_ct1 <- make_scatter_panel(drive, "p_m_join_1", "p_s_join_1", subtitle = "Adipo/Endo", base_font = base_font)
plot_f1_ct2 <- make_scatter_panel(drive, "p_m_join_2", "p_s_join_2", subtitle = "Fibroblast", base_font = base_font)
plot_f1_ct3 <- make_scatter_panel(drive, "p_m_join_3", "p_s_join_3", subtitle = "Epithelial", base_font = base_font)


# Align axes across 4 panels
f1_aligned <- cowplot::align_plots(
  plot_f1_all, plot_f1_ct1, plot_f1_ct2, plot_f1_ct3,
  align = "hv",
  axis = "tblr"
)

figure1 <- cowplot::plot_grid(
  plotlist = f1_aligned,
  ncol = 2,
  labels = "AUTO",
  label_size = 16,
  label_fontface = "bold"
)

print(figure1)

ggsave(
  file.path(figure_dir, "FigureS1_111.pdf"),
  figure1, width = 8, height = 8
)

# ==============================================================================
# 4. Figure 2: QQ (A) + Venn (B)  [keep your existing code below OR paste the
#    Figure 2 section from the previous version unchanged]
# ==============================================================================

# ---- QQ plot (A) from out3 ----
pvals_q    <- out3$p_join
genename_q <- out3$gene_name
mask_q <- is.finite(pvals_q) & !is.na(pvals_q) & !is.na(genename_q) & pvals_q > 0 & pvals_q < 1
p_clean <- pvals_q[mask_q]
g_clean <- as.character(genename_q[mask_q])

p <- pmin(pmax(p_clean, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(p, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)

df_qq <- data.frame(
  obs  = -log10(sort(p_clean)),
  exp  = -log10(ppoints(length(p_clean))),
  gene = g_clean[order(p_clean)]
)
top_genes_df <- head(df_qq, 10)

plot_f2_a <- ggplot(df_qq, aes(x = exp, y = obs)) +
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

# ---- Venn (B): keep geometry + styling; fix singleton patterns ----
n_ct1_only <- length(which(out3$MAP_pattern_nonnull == "100"  & out3$type == "CellTypeSpecific"))
n_ct2_only <- length(which(out3$MAP_pattern_nonnull == "010"  & out3$type == "CellTypeSpecific"))
n_ct3_only <- length(which(out3$MAP_pattern_nonnull == "001"  & out3$type == "CellTypeSpecific"))
n_ct1_ct2  <- length(which(out3$MAP_pattern_nonnull == "110"  & out3$type == "CellTypeSpecific"))
n_ct1_ct3  <- length(which(out3$MAP_pattern_nonnull == "101"   & out3$type == "CellTypeSpecific"))
n_ct2_ct3  <- length(which(out3$MAP_pattern_nonnull == "011"  & out3$type == "CellTypeSpecific"))
n_shared_specific <- length(which(out3$MAP_pattern_nonnull == "111" & out3$type == "CellTypeSpecific"))
n_shared_nonspecific <- length(which(out3$fdr_p_join < 0.1 & out3$type == "NonSpecific"))

ct_labels <- c("Adipo/Endo", "Fibroblast", "Epithelial")
venn_colors <- c("Adipo/Endo" = "#FDB462", "Fibroblast" = "#4E79A7", "Epithelial" = "#59A14F")

circles <- data.frame(
  x0 = c(-1.1,  1.1,  0.0),
  y0 = c( 0.6,  0.6, -1.0),
  r  = c( 2.1,  2.1,  2.1),
  type = ct_labels
)

pos <- list(
  ct1_only = c(-2.25,  0.70),
  ct2_only = c( 2.25,  0.70),
  ct3_only = c( 0.00, -2.40),
  ct1_ct2  = c( 0.00,  1.25),
  ct1_ct3  = c(-0.90, -0.55),
  ct2_ct3  = c( 0.90, -0.55),
  title1   = c(-1.10,  2.80),
  title2   = c( 1.10,  2.80),
  title3   = c( 0.00, -3.20),
  spec_txt = c( 0.00,  0.55),
  line_y   = 0.15,
  nons_txt = c( 0.00, -0.35)
)

plot_f2_b <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = type, color = type),
    alpha = 0.45, linewidth = 0.8
  ) +
  scale_fill_manual(values = venn_colors) +
  scale_color_manual(values = venn_colors) +
  annotate("text", x = pos$ct1_only[1], y = pos$ct1_only[2], label = n_ct1_only,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct2_only[1], y = pos$ct2_only[2], label = n_ct2_only,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct3_only[1], y = pos$ct3_only[2], label = n_ct3_only,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct1_ct2[1], y = pos$ct1_ct2[2], label = n_ct1_ct2,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct1_ct3[1], y = pos$ct1_ct3[2], label = n_ct1_ct3,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct2_ct3[1], y = pos$ct2_ct3[2], label = n_ct2_ct3,
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$title1[1], y = pos$title1[2], label = ct_labels[1],
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$title2[1], y = pos$title2[2], label = ct_labels[2],
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$title3[1], y = pos$title3[2], label = ct_labels[3],
           size = 4, fontface = "bold", color = "black") +
  annotate("text", x = pos$spec_txt[1], y = pos$spec_txt[2],
           label = paste0("Cell Type\nSpecific\n", n_shared_specific),
           size = 2, fontface = "bold", color = "black", lineheight = 0.95) +
  annotate("segment", x = -0.9, xend = 0.9, y = pos$line_y, yend = pos$line_y,
           color = "black", linewidth = 1.2) +
  annotate("text", x = pos$nons_txt[1], y = pos$nons_txt[2],
           label = paste0("Non-specific\n", n_shared_nonspecific),
           size = 2, fontface = "bold", color = "black", lineheight = 0.95) +
  coord_fixed() +
  xlim(-3.5, 3.5) + ylim(-3.2, 3.2) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10))

plot_f2_b_tall <- cowplot::plot_grid(ggdraw(), plot_f2_b, ggdraw(),
                                     ncol = 1, rel_heights = c(0.18, 1, 0.18))

f2_aligned <- cowplot::align_plots(plot_f2_a, plot_f2_b_tall, align = "hv", axis = "tblr")
figure2 <- cowplot::plot_grid(
  f2_aligned[[1]], f2_aligned[[2]],
  ncol = 2, rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 16, label_fontface = "bold",
  label_x = 0, label_y = 1
)

print(figure2)

ggsave(file.path(figure_dir, "FigureS2_111.pdf"),
       figure2, width = 10, height = 5)

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

base_font <- 12

# ==============================================================================
# 1. Load Data
# ==============================================================================
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")
drive  <- read.csv("/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/drive_result_full_lam_new.csv",
                   row.names = 1)

# ==============================================================================
# 2. Plot A: Scatter Plot (classic, no internal "A")
# ==============================================================================
df_c <- drive %>%
  filter(p_m_join > 0, p_m_join < 1, !is.na(p_m_join),
         p_s_join > 0, p_s_join < 1, !is.na(p_s_join)) %>%
  mutate(
    logp_p_m = -log10(p_m_join),
    logp_p_s = -log10(p_s_join)
  )

cor_val   <- cor(df_c$logp_p_m, df_c$logp_p_s, method = "pearson")
cor_label <- paste0("r = ", round(cor_val, 2))
axis_limit <- max(c(df_c$logp_p_m, df_c$logp_p_s), na.rm = TRUE) * 1.05

plot_a <- ggplot(df_c, aes(x = logp_p_m, y = logp_p_s)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.3, size = 1.2) +
  annotate("text", x = 0, y = axis_limit, label = cor_label,
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  xlim(0, axis_limit) + ylim(0, axis_limit) +
  labs(
    x = expression(bold(-log[10](p) ~ "in MiXcan")),
    y = expression(bold(-log[10](p) ~ "in S-MiXcan"))
  ) +
  theme_classic(base_size = base_font) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.line  = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==============================================================================
# 3. Plot B: QQ Plot (classic, no internal "B")
# ==============================================================================
pvals_b    <- merged$p_join_ct3
genename_b <- merged$gene_name
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

plot_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
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
# 5. Plot C: Custom Split Venn (KEEP your hand-written numbers)
# ==============================================================================
# ==============================================================================
# 3 cell types: compute Venn counts (specific via PRIMO MAP + nonspecific count)
# and PLOT a 3-circle Venn using ggforce (classic, no legend).
# ==============================================================================
library(ggplot2)
library(ggforce)
library(cowplot)


# 7 regions
n_ct1_only <- 0
n_ct2_only <- 4
n_ct3_only <- 9
n_ct1_ct2  <- 2
n_ct1_ct3  <- 0
n_ct2_ct3  <- 1
n_shared_specific <- 9

# nonspecific
n_shared_nonspecific <- 8

# -----------------------
# (B) Plot: 3-circle Venn
# -----------------------
label_top    <- "Cell Type\nSpecific"
label_bottom <- "Non-specific"

ct_labels <- c("Cell type 1", "Cell type 2", "Cell type 3")
venn_colors <- c("Cell type 1" = "#FDB462", "Cell type 2" = "#4E79A7", "Cell type 3" = "#59A14F")

# Circle geometry (tweak if you want tighter/looser overlaps)
# ==============================================================================
# FINAL ct3 split-center Venn (ggforce) — tuned to avoid overlaps
# Assumes you already computed:
#   n_ct1_only, n_ct2_only, n_ct3_only,
#   n_ct1_ct2, n_ct1_ct3, n_ct2_ct3, n_ct1_ct2_ct3,
#   n_shared_specific, n_shared_nonspecific
# ==============================================================================

library(ggplot2)
library(ggforce)

# ---- circle labels + colors ----
ct_labels <- c("Adipose", "Fibroblast", "Epithelial")
venn_colors <- c("Adipose" = "#FDB462",
                 "Fibroblast" = "#4E79A7",
                 "Epithelial" = "#59A14F")

# ---- circle geometry (match your screenshot) ----
circles <- data.frame(
  x0 = c(-1.1,  1.1,  0.0),
  y0 = c( 0.6,  0.6, -1.0),
  r  = c( 2.1,  2.1,  2.1),
  type = ct_labels
)

# ---- text positions (hand-tuned) ----
pos <- list(
  ct1_only = c(-2.25,  0.70),
  ct2_only = c( 2.25,  0.70),
  ct3_only = c( 0.00, -2.40),

  ct1_ct2  = c( 0.00,  1.25),
  ct1_ct3  = c(-0.90, -0.55),
  ct2_ct3  = c( 0.90, -0.55),

  triple   = c( 0.00, -0.05),   # slightly below "Cell type 3" label area

  title1   = c(-1.10,  2.80),
  title2   = c( 1.10,  2.80),
  title3   = c( 0.00,  -3.20),   # moved up to avoid overlap

  spec_txt = c( 0.00,  0.55),
  line_y   = 0.15,
  nons_txt = c( 0.00, -0.35)
)

plot_c3_final <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = type, color = type),
    alpha = 0.45, linewidth = 0.8
  ) +
  scale_fill_manual(values = venn_colors) +
  scale_color_manual(values = venn_colors) +

  # ---- region numbers ----
annotate("text", x = pos$ct1_only[1], y = pos$ct1_only[2],
         label = n_ct1_only, size = 5, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct2_only[1], y = pos$ct2_only[2],
           label = n_ct2_only, size = 5, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct3_only[1], y = pos$ct3_only[2],
           label = n_ct3_only, size = 5, fontface = "bold", color = "black") +

  annotate("text", x = pos$ct1_ct2[1], y = pos$ct1_ct2[2],
           label = n_ct1_ct2, size = 4.8, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct1_ct3[1], y = pos$ct1_ct3[2],
           label = n_ct1_ct3, size = 4.8, fontface = "bold", color = "black") +
  annotate("text", x = pos$ct2_ct3[1], y = pos$ct2_ct3[2],
           label = n_ct2_ct3, size = 4.8, fontface = "bold", color = "black") +

  #annotate("text", x = pos$triple[1], y = pos$triple[2],
    #       label = n_ct1_ct2_ct3, size = 4.8, fontface = "bold", color = "black") +

  # ---- titles (all black) ----
annotate("text", x = pos$title1[1], y = pos$title1[2],
         label = ct_labels[1], size = 5, fontface = "bold", color = "black") +
  annotate("text", x = pos$title2[1], y = pos$title2[2],
           label = ct_labels[2], size = 5, fontface = "bold", color = "black") +
  annotate("text", x = pos$title3[1], y = pos$title3[2],
           label = ct_labels[3], size = 5, fontface = "bold", color = "black") +

  # ---- split center ----
annotate("text", x = pos$spec_txt[1], y = pos$spec_txt[2],
         label = paste0("Cell Type\nSpecific\n", n_shared_specific),
         size = 3.8, fontface = "bold", color = "black", lineheight = 0.95) +

  annotate("segment", x = -0.9, xend = 0.9, y = pos$line_y, yend = pos$line_y,
           color = "black", linewidth = 1.2) +

  annotate("text", x = pos$nons_txt[1], y = pos$nons_txt[2],
           label = paste0("Non-specific\n", n_shared_nonspecific),
           size = 3.8, fontface = "bold", color = "black", lineheight = 0.95) +

  coord_fixed() +
  xlim(-3.5, 3.5) +
  ylim(-3.2, 3.2) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

plot_c3_final


# ==============================================================================
# 6. Final Assembly: align A/B axes; make C similar size; then combine
# ==============================================================================

# (1) 强制 A/B 坐标轴区域对齐
ab_aligned <- cowplot::align_plots(
  plot_a, plot_b,
  align = "hv",
  axis = "tblr"
)

plot_ab <- cowplot::plot_grid(
  ab_aligned[[1]], ab_aligned[[2]],
  ncol = 2,
  rel_widths = c(1, 1)
)

# (2) 让 C “视觉高度”接近 A/B（上下加一点空白）
plot_c_tall <- cowplot::ggdraw() +
  cowplot::draw_plot(
    plot_c3_final,
    x = 0, y = 0,
    width = 1, height = 1,
    scale = 1   # 👈 放大倍数，1.15–1.35 都很安全
  )


# (3) 最终拼三块，并统一加 A/B/C label
final_figure <- cowplot::plot_grid(
  ab_aligned[[1]], ab_aligned[[2]], plot_c_tall,
  ncol = 3,
  rel_widths = c(1, 1, 1.2),
  labels = c("A", "B", "C"),
  label_size = 16,
  label_fontface = "bold",
  label_x = 0, label_y = 1,
  align = "hv", axis = "tblr"
)

print(final_figure)

ggsave("/Users/zhusinan/Downloads/FigureS1_5.pdf",
       final_figure, width = 15, height = 5)


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
pvals_b    <- merged$p_join_ct2
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
# 4. Primo Analysis & Filtering (unchanged)
# ==============================================================================

# 2 cell types
# ==============================================================================
# Numbers for Plot C (computed from your provided logic)
# ==============================================================================

# Inputs
fdr_cutoff <- 0.1
specific_label   <- "CellTypeSpecific"
nonspecific_label <- "NonSpecific"

# 1) Split
merged_ctspec <- merged[merged$type_ct2 == specific_label, , drop = FALSE]
merged_unspec <- merged[merged$type_ct2 == nonspecific_label, , drop = FALSE]

# 2) Run PRIMO on ALL specific rows (assume all specific)
res <- primo_map_all(
  merged_ctspec,
  pvals_names = c("p_1_ct2", "p_2_ct2"),
  alt_props   = c(0.05, 0.05)
)

# 3) Significant filter (OR rule), applied within each subset
sig_spec_idx <- which(
  merged_ctspec$fdr_p_1_ct2 < fdr_cutoff | merged_ctspec$fdr_p_2_ct2 < fdr_cutoff
)

sig_uns_idx <- which(
  merged_unspec$fdr_p_1_ct2 < fdr_cutoff | merged_unspec$fdr_p_2_ct2 < fdr_cutoff
)

# 4) MAP class among significant specific rows
#    Using posterior probs from PRIMO result (same row order as merged_ctspec)
post_sig <- res$primo$post_prob[sig_spec_idx, , drop = FALSE]
MAP_class <- max.col(post_sig)

# 5) Assign Plot C numbers (this matches your earlier interpretation)
#    For 2 traits, PRIMO commonly corresponds to patterns:
#    1=null(00), 2=trait1 only(10), 3=trait2 only(01), 4=both(11)
n_trait1 <- sum(MAP_class == 2, na.rm = TRUE)
n_trait2 <- sum(MAP_class == 3, na.rm = TRUE)
n_shared_specific <- sum(MAP_class == 4, na.rm = TRUE)

# 6) Nonspecific significant count (your "n_shared_nonspecific")
n_shared_nonspecific <- length(sig_uns_idx)

# 7) Total shared
n_shared_total <- n_shared_specific + n_shared_nonspecific

# ==============================================================================
# 5. Plot C: Custom Split Venn (KEEP your hand-written numbers)
# ==============================================================================


label_top    <- "Cell Type\nSpecific"
label_bottom <- "Non-specific"

circles <- data.frame(
  x0 = c(-0.8, 0.8),
  y0 = c(0, 0),
  r  = c(2, 2),
  type = c("Epithelial", "Stromal")
)

venn_colors <- c("Epithelial" = "#FDB462", "Stromal" = "#4E79A7")

plot_c_custom <- ggplot() +
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

  # 顶部标题也改成黑色
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
plot_c_tall <- cowplot::plot_grid(
  ggdraw(), plot_c_custom, ggdraw(),
  ncol = 1,
  rel_heights = c(0.18, 1, 0.18)
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

ggsave("/Users/zhusinan/Downloads/Figure2_Final_AJHG_aligned.pdf",
       final_figure, width = 15, height = 5)


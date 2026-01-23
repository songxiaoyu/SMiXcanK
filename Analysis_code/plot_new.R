# ==============================================================================
# 0. Setup: Libraries & Global Settings
# ==============================================================================
library(ggplot2)
library(ggrepel)       # Labels
library(cowplot)       # Final Grid Layout
library(bacon)         # Inflation factor
library(dplyr)
library(tidyr)
library(Primo)
# 专门用于手工画圆的包 (如果没有请安装: install.packages("ggforce"))
library(ggforce)

my_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# 1. Load Data (Please update paths!)
# ==============================================================================
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")
drive <- read.csv('/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/drive_result_full_lam_new.csv', row.names = 1)

# ==============================================================================
# 2. Plot A: Scatter Plot
# ==============================================================================
p_m <- drive$p_m_join
p_s <- drive$p_s_join
df_c <- drive %>%
  filter(p_m > 0, p_m < 1, !is.na(p_m), p_s > 0, p_s < 1, !is.na(p_s)) %>%
  mutate(logp_p_m = -log10(p_m), logp_p_s = -log10(p_s))

cor_val <- cor(df_c$logp_p_m, df_c$logp_p_s, method = "pearson")
cor_label <- paste0("Pearson's r = ", round(cor_val, 2))
axis_limit <- max(c(df_c$logp_p_m, df_c$logp_p_s)) * 1.05

plot_a <- ggplot(df_c, aes(x = logp_p_m, y = logp_p_s)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.3) +
  annotate("text", x = 0, y = axis_limit, label = cor_label, hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  xlim(0, axis_limit) + ylim(0, axis_limit) +
  coord_fixed() +
  labs(title = "A", x = expression(bold(-log[10](p_MiXcan))), y = expression(bold(-log[10](p_SMiXcanK)))) +
  my_theme + theme(panel.grid.major = element_line(color = "grey90"))

# ==============================================================================
# 3. Plot B: QQ Plot
# ==============================================================================
pvals_b <- merged$p_join_ct2
genename_b <- merged$gene_name
mask_b <- is.finite(pvals_b) & !is.na(pvals_b) & !is.na(genename_b) & pvals_b > 0 & pvals_b < 1
p_clean_b <- pvals_b[mask_b]
g_clean_b <- as.character(genename_b[mask_b])

p <- pmin(pmax(p_clean_b, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(1-p, lower.tail = FALSE)
bc <- bacon(y, na.exclude = TRUE)
lambda_val <- inflation(bc)

df_qq <- data.frame(obs = -log10(sort(p_clean_b)), exp = -log10(ppoints(length(p_clean_b))), gene = g_clean_b[order(p_clean_b)])
top_genes_df <- head(df_qq, 8)

plot_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
  geom_text_repel(data = top_genes_df, aes(label = gene), nudge_x = -1, nudge_y = 1, force = 5, box.padding = 0.5, size = 3.5, fontface = "bold") +
  annotate("text", x = 0, y = max(df_qq$obs), label = bquote(lambda[GC] == .(round(lambda_val, 3))), hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  labs(title = "B", x = expression(bold(Expected ~ -log[10](p))), y = expression(bold(Observed ~ -log[10](p)))) +
  my_theme

# ==============================================================================
# 4. Primo Analysis & Filtering
# ==============================================================================
pvals_primo <- merged[, c('p_1_ct2', 'p_2_ct2')]
colnames(pvals_primo) <- c("P1", "P2")
rownames(pvals_primo) <- merged$gene_name

res <- Primo::Primo_pval(pvals = pvals_primo, alt_props = c(0.05, 0.05))

fdr_trait1 <- p.adjust(pvals_primo$P1, method = "BH")
fdr_trait2 <- p.adjust(pvals_primo$P2, method = "BH")
map_index <- max.col(res$post_prob)
cutoff <- 0.1

# 使用 OR 逻辑筛选 (只要任意一侧显著就算)
final_trait1_genes <- rownames(pvals_primo)[which(map_index == 2 & (fdr_trait1 < cutoff | fdr_trait2 < cutoff))]
final_trait2_genes <- rownames(pvals_primo)[which(map_index == 3 & (fdr_trait1 < cutoff | fdr_trait2 < cutoff))]
final_shared_genes <- rownames(pvals_primo)[which(map_index == 4 & (fdr_trait1 < cutoff | fdr_trait2 < cutoff))]

# ==============================================================================
# 5. Plot C: Custom Split Venn (Ultimate Version)
# ==============================================================================
# 基础计数
n_trait1 <- length(final_trait1_genes)
n_trait2 <- length(final_trait2_genes)
n_shared_total <- length(final_shared_genes)
n_trait1 <- 4
n_trait2 <- 8
n_shared_total <- 64


# [示例逻辑]：请替换为你真实的计算逻辑
n_shared_specific    <- 59
n_shared_nonspecific <- 5

# 定义标签名称
label_top    <- "Cell Type\nSpecific"
label_bottom <- "Non-specific"
# -----------------------------------

# 准备画图数据
circles <- data.frame(
  x0 = c(-0.8, 0.8),
  y0 = c(0, 0),
  r  = c(2, 2),
  type = c("cell type 1", "cell type 2")
)

# 颜色映射
venn_colors <- c("cell type 1" = "#FDB462", "cell type 2" = "#4E79A7")

# 计算总 N (用于百分比)
N_total <- n_trait1 + n_trait2 + n_shared_total

plot_c_custom <- ggplot() +
  # 1. 画圆 (ggforce)
  geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = type, color = type),
              alpha = 0.5, linewidth = 0.5) +
  scale_fill_manual(values = venn_colors) +
  scale_color_manual(values = venn_colors) + # 边框同色，或者设为 "white"

  # 2. 画交集横线 (白色分割线)
  geom_segment(aes(x = -0.95, xend = 0.95, y = 0, yend = 0),
               color = "white", linewidth = 1.2) +

  # 3. 标注文字
  # 左 (Trait 1 Only)
  annotate("text", x = -1.8, y = 0,
           label = paste0(n_trait1),
           size = 5, fontface = "bold") +

  # 右 (Trait 2 Only)
  annotate("text", x = 1.8, y = 0,
           label = paste0(n_trait2),
           size = 5, fontface = "bold") +

  # 中上 (Specific)
  annotate("text", x = 0, y = 0.7,
           label = paste0(label_top, "\n", n_shared_specific),
           size = 3.5, fontface = "bold", color = "white", lineheight=0.9) +

  # 中下 (Non-specific)
  annotate("text", x = 0, y = -0.7,
           label = paste0(label_bottom, "\n", n_shared_nonspecific),
           size = 3.5, fontface = "bold", color = "white", lineheight=0.9) +

  # 顶部标题 (Cell Type 1 / 2)
  annotate("text", x = -0.8, y = 2.4, label = "cell type 1", size = 5, fontface = "bold", color = "#FDB462") +
  annotate("text", x = 0.8, y = 2.4, label = "cell type 2", size = 5, fontface = "bold", color = "#4E79A7") +

  # 4. 样式调整
  coord_fixed() +
  labs(title = "C") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0, margin = margin(b = 10)),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# ==============================================================================
# 6. Final Assembly (A B C Side-by-Side)
# ==============================================================================
final_figure <- plot_grid(
  plot_a, plot_b, plot_c_custom,
  ncol = 3, nrow = 1,
  # 给 C 图稍微多一点宽度，因为它比较宽
  rel_widths = c(1, 1, 1.2),
  align = 'h', axis = 'tb'
)

# 添加总标题
final_figure_with_title <- plot_grid(
  ggdraw(),
  final_figure,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

print(final_figure_with_title)

# 保存文件
ggsave("/Users/zhusinan/Downloads/Figure2_Final_SplitShared.pdf", final_figure_with_title, width = 17, height = 6)


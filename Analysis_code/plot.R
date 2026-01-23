# ==============================================================================
# Prerequisites: Install necessary packages if missing
# install.packages(c("ggplot2", "ggrepel", "ggVennDiagram", "cowplot", "bacon", "dplyr", "tidyr"))
# ==============================================================================

library(ggplot2)
library(ggrepel)       # For non-overlapping labels in QQ plot
library(ggVennDiagram) # For Venn Diagram
library(cowplot)       # For arranging the final grid
library(bacon)         # For genomic inflation estimation
library(dplyr)         # For data manipulation
library(tidyr)

# Set a consistent theme for all plots
my_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# SECTION 0: Load Your Data
# ==============================================================================
merged <- read.csv("/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/merged_ct3_ct2.csv")
drive <- read.csv('/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/drive_result_full_lam_new.csv', row.names = 1)

# ==============================================================================
# PLOT A: P-value Comparison Scatter Plot (MiXcan & SMiXcan)
# ==============================================================================

# --- 1. Data Extraction & Transformation ---
# Filter for valid p-values in both columns to avoid log(0) or log(NA)
p_m <- drive$p_m_join
p_s <- drive$p_s_join
df_c <- drive %>%
  filter(p_m  > 0, p_m < 1, !is.na(p_m),
         p_s > 0, p_s < 1, !is.na(p_s)) %>%
  mutate(
    logp_p_m = -log10(p_m),
    logp_p_s = -log10(p_s)
  )

# --- 2. Calculate Correlation ---
cor_val <- cor(df_c$logp_p_m, df_c$logp_p_s, method = "pearson")
cor_label <- paste0("Pearson's r = ", round(cor_val, 2))

# Determine axis limits for a square plot
axis_limit <- max(c(df_c$logp_p_m, df_c$logp_p_s)) * 1.05

# --- 3. Plotting ---
plot_a <- ggplot(df_c, aes(x = logp_p_m, y = logp_p_s)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_point(color = "#4E79A7", alpha = 0.3) +
  annotate("text", x = 0, y = axis_limit, label = cor_label,
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  xlim(0, axis_limit) + ylim(0, axis_limit) +
  coord_fixed() +
  labs(title = "A",
       x = expression(bold(-log[10](p_MiXcan))),
       y = expression(bold(-log[10](p_SMiXcanK)))) +
  my_theme + theme(panel.grid.major = element_line(color = "grey90"))


# ==============================================================================
# PLOT B: S-MiXcan_K QQ Plot (Using p_join_ct2)
# ==============================================================================

# --- 1. Data Extraction & Cleaning ---
pvals_b <- merged$p_join_ct2
genename_b <- merged$gene_name

# Clean data: ensure finite, non-NA, and strictly within (0,1) for qnorm
mask_b <- is.finite(pvals_b) & !is.na(pvals_b) & !is.na(genename_b) & pvals_b > 0 & pvals_b < 1
p_clean_b <- pvals_b[mask_b]
g_clean_b <- as.character(genename_b[mask_b])

# --- 2. Bacon Calculation for Genomic Inflation ---
p <- pmin(pmax(p_clean_b, .Machine$double.eps), 1 - .Machine$double.eps)
y <- qnorm(1-p, lower.tail = FALSE)

bc <- bacon(y, na.exclude = TRUE)
estimates(bc)
lambda_val <- inflation(bc)

# --- 3. Prepare Plotting Data Frame ---
df_qq <- data.frame(
  obs = -log10(sort(p_clean_b)),
  exp = -log10(ppoints(length(p_clean_b))),
  gene = g_clean_b[order(p_clean_b)]
)

# Select top genes to label (e.g., top 3)
top_genes_df <- head(df_qq, 8)

# --- 4. Plotting ---
# plot_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
#   geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +
#   # Label top genes with ggrepel
#   geom_text_repel(data = top_genes_df, aes(label = gene),
#                   nudge_x = -1, nudge_y = 1, force = 5,
#                   box.padding = 0.5, point.padding = 0.2,
#                   arrow = arrow(length = unit(0.02, "npc")),
#                   size = 3.5, fontface = "bold") +
#   # Annotation for Lambda
#   annotate("text", x = max(df_qq$exp), y = 0,
#            label = bquote(lambda[GC] == .(round(lambda_val, 3))),
#            hjust = 1, vjust = 0, size = 5, fontface = "bold") +
#   labs(title = "B",
#        x = expression(bold(Expected ~ -log[10](p))),
#        y = expression(bold(Observed ~ -log[10](p)))) +
#   my_theme

plot_b <- ggplot(df_qq, aes(x = exp, y = obs)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  geom_point(color = "#4E79A7", alpha = 0.8, size = 2) +

  # Label top genes with ggrepel
  geom_text_repel(data = top_genes_df, aes(label = gene),
                  nudge_x = -1, nudge_y = 1, force = 5,
                  box.padding = 0.5, point.padding = 0.2,
                  arrow = arrow(length = unit(0.02, "npc")),
                  size = 3.5, fontface = "bold") +

  # --- 修改开始: 将 Lambda 注释移到左上角 ---
  annotate("text",
           x = 0,                  # x设为0 (最左侧)
           y = max(df_qq$obs),     # y设为最大观测值 (最顶部)
           label = bquote(lambda[GC] == .(round(lambda_val, 3))),
           hjust = 0,              # 0表示左对齐 (Left alignment)
           vjust = 1,              # 1表示顶对齐 (Top alignment)
           size = 5, fontface = "bold") +
  # --- 修改结束 ---

  labs(title = "B",
       x = expression(bold(Expected ~ -log[10](p))),
       y = expression(bold(Observed ~ -log[10](p)))) +
  my_theme

# ==============================================================================
# PLOT A: Regulatory Pattern Proportion (Stacked Bar Chart)
# Based on your Primo analysis code
# ==============================================================================

# --- 1. Run Your Primo Analysis (You need to run this part first) ---
# Ensure 'combined2' and 'Primo' package are loaded before this step.
# ... your code to load data ...
# library('Primo')
# pvals2 = merged[,c('p_1_ct2','p_2_ct2')]
# primo_results2 <- Primo_pval(pvals=pvals2, alt_props=c(0.05, 0.05))
# pp <- primo_results2$post_prob
# rownames(pp) <- merged$gene_name
# colnames(pp) <- c('00','10', '01', '11')
#
# # --- 2. Process Primo Results into a Data Frame for Plotting ---
# # 1) assign each gene to its highest-probability pattern
# best_idx <- max.col(pp, ties.method = "first")
# best_pat <- colnames(pp)[best_idx]
#
# # 2) counts + proportions (keep all patterns, even if 0)
# pat_levels <- colnames(pp)
# # Use factor with levels to ensure all patterns are present even if count is 0
# tab <- table(factor(best_pat, levels = pat_levels))
# prop <- as.numeric(tab) / sum(tab)
#
# # Create the data frame for plotting
# df_a <- data.frame(
#   pattern = pat_levels,
#   prop = prop
# )
#
# # (Optional) Drop '00' pattern and re-normalize proportions
# # This focuses the chart on the composition of non-null associations
# df_a <- subset(df_a, pattern != "00")
# df_a$prop <- df_a$prop / sum(df_a$prop)
#
# # Calculate cumulative proportions and midpoints for label placement
# df_a <- df_a %>%
#   arrange(desc(pattern)) %>%
#   mutate(
#     cumulative_prop = cumsum(prop),
#     midpoint = cumulative_prop - (prop / 2),
#     # Create a clear label with pattern and percentage
#     label_text = paste0(pattern, "\n(", round(prop * 100, 1), "%)")
#   )
#
# # --- 3. Plotting with ggplot2 ---
#
# # Define a color palette that is distinct and visually appealing
# # 11=Purple, 10=Green, 01=Orange (Example colors, adjust as you like)
# my_colors_a <- c("11" = "#984EA3", "10" = "#4DAF4A", "01" = "#FF7F00")
#
# plot_c <- ggplot(df_a, aes(x = 1, y = prop, fill = pattern)) +
#   # Create the stacked bar
#   geom_col(width = 0.7, color = "white", size = 0.2) +
#   # Add text labels at the midpoint of each segment
#   geom_text(aes(y = midpoint, label = label_text),
#             color = "white", size = 3.5, fontface = "bold") +
#   # Make the plot horizontal
#   coord_flip() +
#   # Apply custom colors
#   scale_fill_manual(values = my_colors_a) +
#   # Remove x-axis ticks and labels (since it's just one bar)
#   scale_x_continuous(breaks = NULL) +
#   # Format the y-axis as percentages
#   scale_y_continuous(labels = scales::percent_format()) +
#   # Add titles and labels
#   labs(title = "C",
#        x = NULL,
#        y = "Proportion of Non-Null Genes",
#        fill = "Association Pattern") +
#   # Use a clean theme
#   theme_minimal(base_size = 12) +
#   theme(
#     plot.title = element_text(face = "bold"),
#     panel.grid.major.y = element_blank(), # Remove horizontal grid lines
#     panel.grid.minor = element_blank(),
#     axis.text.y = element_blank() # Remove y-axis text (the '1')
#   )
counts_data <- data.frame(
  pattern = c("10", "01", "11"),
  count = c(length(final_trait1_genes),  # Trait 1 Specific (e.g., 6)
            length(final_trait2_genes),  # Trait 2 Specific (e.g., 11)
            length(final_shared_genes))  # Shared (e.g., 40)
)

# 计算比例和标签位置
df_plot <- counts_data %>%
  # 只有当总数大于0时才计算，防止报错
  filter(count > 0) %>%
  mutate(
    prop = count / sum(count),
    # 重新排序 pattern 以便画图顺序好看 (10 -> 11 -> 01)
    pattern = factor(pattern, levels = c("10", "11", "01")),
    # 计算累积比例用于定位标签
    cumulative_prop = cumsum(prop),
    midpoint = cumulative_prop - (prop / 2),
    # 标签内容：类别 + 数量 + 百分比
    label_text = paste0(pattern, "\nN=", count, "\n(", percent(prop, accuracy = 0.1), ")")
  )

# ==============================================================================
# 2. 设置颜色 (为了配合之前的韦恩图和 QQ 图)
# ==============================================================================
# 10 (Trait 1) = 红色系
# 01 (Trait 2) = 蓝色系
# 11 (Shared)  = 紫色系
my_colors_c <- c("10" = "#FDB462",  # 对应 Trait 1 (接近橙/红)
                 "01" = "#4E79A7",  # 对应 Trait 2 (蓝色)
                 "11" = "#984EA3")  # Shared (紫色)

# ==============================================================================
# 3. 绘图 (Figure C)
# ==============================================================================
plot_c <- ggplot(df_plot, aes(x = 1, y = prop, fill = pattern)) +
  # 1. 画堆叠柱状图
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +

  # 2. 添加文字标签 (N 和 %)
  geom_text(aes(y = midpoint, label = label_text),
            color = "white", size = 4, fontface = "bold") +

  # 3. 横向翻转 (变成条形图)
  coord_flip() +

  # 4. 颜色设置
  scale_fill_manual(values = my_colors_c,
                    labels = c("10" = "Trait 1 Specific",
                               "01" = "Trait 2 Specific",
                               "11" = "Shared")) +

  # 5. 去除多余的轴刻度
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(labels = percent_format()) +

  # 6. 标题和主题
  labs(title = "C",
       subtitle = paste0("Distribution of Significant Genes (Total N=", sum(counts_data$count), ")"),
       x = NULL,
       y = "Proportion",
       fill = "Association Pattern") +

  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),     # 去掉所有网格线
    axis.text.y = element_blank(),    # 去掉Y轴文字
    legend.position = "bottom"        # 图例放下面
  )
# ==============================================================================
# PLOT D: Significant Gene Overlap (p_join_ct2 vs p_join_ct3)
# ==============================================================================

# --- 1. Define Threshold and Extract Gene Lists ---

fdr_trait1 <- p.adjust(merged_1$p_1_ct2, method = "BH")
fdr_trait2 <- p.adjust(merged_1$p_2_ct2, method = "BH")
fdr_shared <- p.adjust(merged_1$p_join_ct2, method = "BH")

# 2. 设定标准
cutoff <- 0.1

# 3. 最终筛选 (Intersection)
# 必须满足：Primo 说是特异性 AND 统计学显著
map_index <- max.col(res$post_prob)
final_trait1_genes <- rownames(pvals)[which(map_index == 2 & fdr_trait1 < cutoff)]
final_trait2_genes <- rownames(pvals)[which(map_index == 3 & fdr_trait2 < cutoff)]
final_shared_genes <- rownames(pvals)[which(map_index == 4 & fdr_trait1 < cutoff & fdr_trait2 < cutoff)]

merged$fwer_p_ct2
# 查看数量
length(final_trait1_genes)
length(final_trait2_genes)
length(final_shared_genes)
# genes_ct3 <- merged$gene_name[which(merged$fwer_p_ct3 < sig_thresh_venn)]
# genes_ct2 <- merged$gene_name[which(merged$fwer_p_ct2 < sig_thresh_venn)]
#
# # --- 2. Prepare Data List for ggVennDiagram ---
# venn_data <- list(
#   "cell type 1" = final_trait1_genes,
#   "cell type 2" = final_trait2_genes
# )
#
# # --- 3. Plotting ---
# # Using a blue-to-orange gradient to match the theme
# plot_d <- ggVennDiagram(venn_data,
#                         label_alpha = 0,
#                         category.names = c("cell type 1", "cell type 2"),
#                         set_size = 5, edge_size = 0) +
#   scale_fill_gradient(low = "#FDB462", high = "#4E79A7") +
#   labs(title = paste0("D")) +
#   theme_void() +
#   theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
#         legend.position = "none")

list_trait1 <- unique(c(final_trait1_genes, final_shared_genes))

# 集合 B: Trait 2 的所有显著基因 (特异 + 共享)
list_trait2 <- unique(c(final_trait2_genes, final_shared_genes))

# --- 2. 准备 ggVennDiagram 的输入列表 ---
venn_data <- list(
  "Trait 1" = list_trait1,
  "Trait 2" = list_trait2
)

# --- 3. 画图 ---
plot_d <- ggVennDiagram(venn_data,
                        label_alpha = 0,
                        category.names = c("Trait 1", "Trait 2"),
                        set_size = 5,
                        edge_size = 0) +
  # 设置颜色渐变 (蓝到橙)
  scale_fill_gradient(low = "#FDB462", high = "#4E79A7") +
  labs(title = "D. Shared & Specific Genes") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "none")

# ==============================================================================
# FINAL ASSEMBLY: Combine into a 2x2 Grid
# ==============================================================================

final_figure <- plot_grid(
  plot_a, plot_b, plot_c, plot_d,
  ncol = 2, nrow = 2,
  align = 'hv', axis = 'tblr'
)

# Add main title
final_figure_with_title <- plot_grid(
  ggdraw() + draw_label("Figure 2", fontface='bold', size=16),
  final_figure,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Display plots
print(final_figure_with_title)

# Save output (optional)
ggsave("/Users/zhusinan/Downloads/Figure2_with_merged_data4.pdf", final_figure_with_title, width = 12, height = 10)


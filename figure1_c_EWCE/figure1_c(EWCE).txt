ibrary(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix)
###Step 1: 准备数据###
##scRNA_combined为整合单细胞数据
##计算Lvl4平均表达量
avg_exp <- AverageExpression(scRNA_combined, 
                             assays = "RNA",
                             group.by = "cells",
                             return.seurat = TRUE)

# 提取标准化表达矩阵
exp_matrix <- Seurat::GetAssayData(
  avg_exp,
  assay = "RNA",
  layer = "counts"
)

# 转换为长格式
exp_lvl5 <- exp_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Lvl4",
    values_to = "Expr_sum_mean"
  )
##过滤未表达基因
not_expressed <- exp_lvl5 %>%
  group_by(Gene) %>%
  summarise(total_sum = sum(Expr_sum_mean)) %>%
  filter(total_sum == 0) %>%
  pull(Gene)

exp_lvl5 <- exp_lvl5 %>%
  filter(!Gene %in% not_expressed)

##按细胞类型标准化
exp_lvl5 <- exp_lvl5 %>%
  group_by(Lvl4) %>%
  mutate(Expr_sum_mean = Expr_sum_mean * 1e6 / sum(Expr_sum_mean))

##计算基因特异性
exp_lvl5 <- exp_lvl5 %>%
  group_by(Gene) %>%
  mutate(specificity = Expr_sum_mean / sum(Expr_sum_mean)) %>%
  ungroup()

#exp_lvl5 <- readRDS("/data/h01005/hc/ldsc/exp_lvl5.RDS")

ewce_analysis <- function(target_genes, 
                          background_genes, 
                          specificity_df, 
                          n_perm = 1000,
                          p_adjust_method = "BH") { # 默认改为FDR控制
  
  # 列名标准化
  colnames(specificity_df) <- c("Gene", "CellType", "Expr_sum_mean", "specificity")
  
  # 过滤基因
  target_genes <- intersect(target_genes, specificity_df$Gene)
  background_genes <- intersect(background_genes, specificity_df$Gene)
  
  # 计算实际观测值（关键修正）
  obs_scores <- specificity_df %>%
    filter(Gene %in% target_genes) %>%
    group_by(CellType) %>%
    summarise(obs_score = mean(specificity), .groups = "drop")
  
  # 初始化矩阵
  cell_types <- unique(specificity_df$CellType)
  perm_matrix <- matrix(nrow = n_perm, ncol = length(cell_types),
                        dimnames = list(NULL, cell_types))
  
  # 置换检验
  for(i in 1:n_perm){
    random_genes <- sample(background_genes, length(target_genes))
    temp_scores <- specificity_df %>%
      filter(Gene %in% random_genes) %>%
      group_by(CellType) %>%
      summarise(score = mean(specificity), .groups = "drop")
    
    # 确保细胞类型顺序一致
    perm_matrix[i, ] <- temp_scores$score[match(cell_types, temp_scores$CellType)]
  }
  
  # 计算p值（向量化计算）
  results <- obs_scores %>%
    mutate(
      p_value = sapply(1:nrow(.), function(j) {
        ct <- .$CellType[j]
        mean(perm_matrix[, ct] >= .$obs_score[j])
      }),
      p_adj = p.adjust(p_value, method = p_adjust_method)
    )
  
  return(results)
}

### Step 4: 执行分析（需先转换数据格式）###
# 定义目标基因和背景基因
target_genes <- c("OPALIN", "MOG", "MBP", "MOBP")
# 执行两组间差异分析
Idents(scRNA_combined) <- "group"
markers_after <- FindMarkers(
  object = scRNA_combined,
  ident.1 = "pd",    # 对照组
  ident.2 = "con",     # 实验组
  only.pos = FALSE,   # 获取双向差异基因
  min.pct = 0.25,     # 基因在至少25%细胞中表达
  logfc.threshold = 0.25,
  test.use = "wilcox",# 默认Wilcoxon秩和检验
  verbose = TRUE
)
# 添加基因名列
markers_after$gene <- rownames(markers_after)
#saveRDS(markers_after,"/data/h01005/hc/再整合-0.1ldsc-多gwas数据/EWCE/markers_after(全细胞).rds")
#提取上调，下调基因
target_genes_up <- markers_after %>%
  slice_max(order_by = avg_log2FC, n = 600)
target_genes_up <- na.omit(target_genes_up)
target_genes_up <- target_genes_up$gene

target_genes_down <- markers_after %>%
  slice_min(order_by = avg_log2FC, n = 600)
target_genes_down <- na.omit(target_genes_down)
target_genes_down <- target_genes_down$gene

background_genes <- unique(exp_lvl5$Gene)  # 使用数据中所有基因作为背景

# 运行分析
results_up <- ewce_analysis(
  target_genes = target_genes_up,
  background_genes = background_genes,
  specificity_df = exp_lvl5
)

results_down <- ewce_analysis(
  target_genes = target_genes_down,
  background_genes = background_genes,
  specificity_df = exp_lvl5
)

### Step 5: 优化后的可视化 ###

###改进##
# 合并结果并添加方向标识
results_combined <- rbind(
  results_up %>% mutate(Direction = "Upregulated", obs_score = obs_score),
  results_down %>% mutate(Direction = "Downregulated", obs_score = -obs_score)
)


# 创建细胞类型顺序映射（保持上下调结果对齐）
cell_order <- results_combined %>%
  arrange(desc(obs_score)) %>%
  distinct(CellType) %>%
  pull(CellType)

# 转换数据格式
plot_data <- results_combined %>%
  mutate(
    CellType = factor(CellType, levels = rev(cell_order)),
    Direction = factor(Direction, levels = c("Upregulated", "Downregulated")),
    sig_label = ifelse(p_adj < 0.05, "*", "")
  )

plot_data <- results_combined %>%
  mutate(
    CellType = factor(CellType, levels = cell_order),
    Direction = factor(Direction, levels = c("Upregulated", "Downregulated")),
    sig_label = ifelse(p_adj < 0.05, "*", ""),
    color_group = case_when(
      p_adj < 0.05 & Direction == "Upregulated" ~ "up_sig",
      p_adj < 0.05 & Direction == "Downregulated" ~ "down_sig",
      TRUE ~ "ns"
    )
  )


# 绘图代码
ggplot(plot_data, aes(y = obs_score, x = CellType, fill = color_group)) +
  geom_col(width = 0.7) +
  # 添加分割线（在y=0位置）
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "dashed") +
  geom_text(
    aes(label = sig_label,
        # 调整星号位置逻辑
        y = ifelse(Direction == "Upregulated", 
                   obs_score + max(abs(plot_data$obs_score))*0.02,  # 上调方向向上偏移
                   obs_score - max(abs(plot_data$obs_score))*0.02)), # 下调方向向下偏移
    color = "black",
    size = 5,
    vjust = ifelse(plot_data$Direction == "Upregulated", 0, 1) # 调整垂直对齐
  ) +
  scale_y_continuous(
    name = "Mean Specificity Score",
    breaks = seq(-1, 1, 0.5),
    labels = abs(seq(-1, 1, 0.5)),
    limits = c(-1, 1) * max(abs(plot_data$obs_score)) * 1.1
  ) +
  scale_fill_manual(
    values = c("up_sig" = "#E41A1C", "down_sig" = "#377EB8", "ns" = "grey80"),
    breaks = c("up_sig", "down_sig", "ns"),
    labels = c("Upregulated", "Downregulated", "Not significant")
  ) +
  labs(title = "Cell Type Enrichment Analysis",
       subtitle = "Upregulated vs Downregulated genes",
       fill = "Significance") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

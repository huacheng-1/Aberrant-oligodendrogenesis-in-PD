library(Seurat)
library(AUCell)
library(GSEABase)
library(dplyr)
library(ggplot2)
library(tidydr)
library(gplots)
library(tidyverse)
library(cowplot)
library(ggthemes)
library(clusterProfiler)

##读取整合数据，生成cells_rankings分数
scRNA_combined <- readRDS("/data/h01005/hc/再整合/scRNA_combined.RDS")
scRNA_pdcon <- subset(scRNA_combined, subset = group %in% c("pd","con"))
rm(scRNA_combined)
#saveRDS(scRNA_pdcon, "/data/h01005/hc/AUCell/20251202/scRNA_pdcon.RDS")
gc()

metadata <- scRNA_pdcon@meta.data
#saveRDS(metadata, "/data/h01005/hc/AUCell/20251202/metadata_pdcon.RDS")

# 转换为普通矩阵（AUCell需要）
expr_matrix <- scRNA_pdcon@assays$SCT@data  # 默认使用"RNA" assay的log-normalized数据
expr_matrix <- as.matrix(expr_matrix)
#saveRDS(expr_matrix, "/data/h01005/hc/AUCell/20251202/expr_matrix_pdcon.RDS")

gc()
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = FALSE)
#saveRDS(cells_rankings, "/data/h01005/hc/AUCell/20251202/cells_rankings_pdcon.RDS")
gc()

#计算基因排名
##这是两个数据集的目标基因集，用于不同数据集的分析。因为处理步骤相同故放在一起

#GSE7621基因（氧化压力）
oxidative_genes_vector1 <- c("ATF4","CAPN2","BANF1","PNPT1","PTPRK","FOXO3","SIRPA","STAU1","ETV5","RELA","KEAP1","ATM","DHFR","ABL1")
#GSE49036基因(氧化压力)
oxidative_genes_vector2 <- c("NFE2L2","CYP1B1","HMOX1","PTGS2","TXNIP","ATM","CRYAB","PTPRK","KAT2B","DDR2","CAT","GPR37" )

oxidative_genes_vector <- list(oxidative_genes_vector1 =oxidative_genes_vector1,
                               oxidative_genes_vector2 =oxidative_genes_vector2)
#画图
for (i in names(oxidative_genes_vector)){
  print(i)
  cells_AUC <- AUCell_calcAUC(oxidative_genes_vector[[i]], cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
  auc_matrix <- getAUC(cells_AUC) 
  head(auc_matrix[, 1:5])
  auc_df <- data.frame(
    Cell = colnames(auc_matrix),
    AUC = auc_matrix[1, ],
    Cluster = metadata$cells,  # 假设metadata中有聚类信息
    Group = metadata$group
  )
  #auc_con <- subset(auc_df, Group == "con")
  #auc_pd <- subset(auc_df, Group == "pd")
  #箱线图
  # 计算每对比较的效应量（如均值差异）
  auc_summary <- auc_df %>%
    group_by(Cluster) %>%
    summarise(mean_AUC = mean(AUC))
  #按均值排序簇
  cluster_order <- auc_summary %>% arrange(mean_AUC) %>% pull(Cluster)
  auc_df$Cluster <- factor(auc_df$Cluster, levels = cluster_order)
  write.csv(auc_df, paste0("/data/h01005/hc/AUCell/20251221/",i,".csv"))
  p <- 
    ggplot(auc_df, aes(x = Cluster, y = AUC, fill = Cluster)) +
    geom_boxplot(outlier.size = 0.5) +
    labs(title = paste0("AUC Scores Across Cell Clusters in PD"), y = "AUC Value") +
    scale_fill_manual(values = c("Oligodendrocyte" = "red", "Astrocyte" = "#F78000", "Oligodendrocyte precursor cell" = "yellow",
                                 "Neuron" ="pink" ,"Microglial cell" = "#3176B7","Fibroblast" = "brown")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  saveRDS(p, paste0("/data/h01005/hc/AUCell/20251221/",i,".RDS"))
  }


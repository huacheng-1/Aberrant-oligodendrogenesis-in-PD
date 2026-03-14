####################第六步-在r中计算p值##############################
library("tidyverse")
library("stringr")
files <- list.files("/data/h01005/hc/ldsc/ldsc/",pattern=".results",full.names = TRUE)
#files <- files[grepl("age",files)]

d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
         makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest() %>% 
  filter(Category=="L2_0") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>% 
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  dplyr::select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)

write_tsv(d,path="/data/h01005/hc/ldsc/ldsc/cell_types_GWAS_pvalues.txt")
#################################################################

#################第七步-p值画图（figure1_a&b）########################
library(tidyr)  
library(ggplot2)
library(tidyverse)
###scRNA###
finall_result <- data.frame(
  celltype = c("astrocyte", "fibroblast", "microglial_cell", "neuron", "oligodendrocyte", "oligodendrocyte_precursor_cell"),
  ieu_b_7 = c(0.46841380097643104, 0.4599130283932974, 0.5849195730947503, 0.5041866058564755, 0.36373987360744664, 0.7110644408746453),
  LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 = c(0.5454709150587502, 0.8365303856443889, 0.022398566647692908, 0.45625740013938076, 0.004665446797138428, 0.2837346315062562),
  GP2_EUR_ONLY = c(0.6485323235785059, 0.9755694009841253, 0.03600845652966844, 0.3378194085179831, 0.0350128676757504, 0.1752248710979334),
  Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry = c(0.2019960321583696, 0.7289435573082299, 0.04621304519585989, 0.7803115727872117, 0.13805993242476333, 0.1016093796315225) 
)
finall_result$ieu_b_7 <- -log10(finall_result$ieu_b_7)
finall_result$LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 <- -log10(finall_result$LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38)
finall_result$GP2_EUR_ONLY <- -log10(finall_result$GP2_EUR_ONLY)
finall_result$Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry <- -log10(finall_result$Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry)

# 转换数据为长格式以便绘图  
library(tidyr)
library(ggplot2)
plot_data <- pivot_longer(finall_result,   
                          cols = c("ieu_b_7", "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
                                   "GP2_EUR_ONLY","Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"),  
                          names_to = "Dataset",  
                          values_to = "Value")  

# 1. 数据预处理：确保组合唯一性并显式排序
plot_data_clean <- plot_data %>%
  # 去除重复组合（保留最大值）
  group_by(celltype, Dataset) %>%
  summarise(Value = max(Value), .groups = "drop") %>%
  # 显式定义因子顺序（先按细胞类型，再按数据集）
  mutate(
    celltype = factor(celltype, 
                      levels = unique(celltype[order(celltype)]),
                      ordered = TRUE),
    Dataset = factor(Dataset,
                     levels = c("ieu_b_7", "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
                                "GP2_EUR_ONLY","Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"),
                     ordered = TRUE)
  ) %>%
  # 物理排序数据框（确保因子顺序与数据行顺序一致）
  arrange(celltype, Dataset)

ggplot(plot_data_clean, aes(x = celltype, y = Value, fill = Dataset)) +  
  geom_col(position = "dodge", width = 0.7) +  # 使用dodge位置并调整柱宽  
  geom_text(aes(label = round(Value, 3)),   
            position = position_dodge(width = 0.7),  
            vjust = -0.5, size = 3) +  
  scale_fill_manual(values = c("#1f78b4", "#33a02c","pink","#f57c6e"),  
                    labels = c("ieu_b_7", "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
                               "GP2_EUR_ONLY","Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry")) +  
  labs(title = "sc_RNA Cell Type Comparison Across GWAS Datasets",  
       x = "Cell Type",  
       y = "-log10(pval)") +  
  theme_minimal() +  
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    legend.title = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  
  ) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 调整y轴范围  


######ATAC####
library(tidyr)  
library(ggplot2)
library(tidyverse)
finall_result <- data.frame(
  celltype = c("astrocyte", "microglial_cell", "neuron", "oligodendrocyte", "oligodendrocyte_precursor_cell"),
  ieu_b_7 = c(0.2069535199204191, 0.35434086362761463, 0.2803916920577668, 0.39677006108164237, 0.4822421736131449),
  LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 = c(0.09134056287402958, 0.80436976630058, 0.6718499714043197, 0.0190978322629749, 0.18286522849836762),
  GP2_EUR_ONLY = c(0.8329032069409021, 0.6563379317056216, 0.032279641914275126, 0.6979193239835468, 0.6571429023425331),
  Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry = c(0.15446340636008493, 0.02110552006675248, 0.01636324551245405, 0.07445546398896186, 0.061631625319569894)
)
finall_result$ieu_b_7 <- -log10(finall_result$ieu_b_7)
finall_result$LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38 <- -log10(finall_result$LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38)
finall_result$GP2_EUR_ONLY <- -log10(finall_result$GP2_EUR_ONLY)
finall_result$Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry <- -log10(finall_result$Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry)
# 转换数据为长格式以便绘图  
plot_data <- pivot_longer(finall_result,   
                          cols = c("ieu_b_7", "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38","GP2_EUR_ONLY", "Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"),  
                          names_to = "Dataset",  
                          values_to = "Value")  


# 1. 数据预处理：确保组合唯一性并显式排序
plot_data_clean <- plot_data %>%
  # 去除重复组合（保留最大值）
  group_by(celltype, Dataset) %>%
  summarise(Value = max(Value), .groups = "drop") %>%
  # 显式定义因子顺序（先按细胞类型，再按数据集）
  mutate(
    celltype = factor(celltype, 
                      levels = unique(celltype[order(celltype)]),
                      ordered = TRUE),
    Dataset = factor(Dataset,
                     levels = c("ieu_b_7", 
                                "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
                                "GP2_EUR_ONLY",
                                "Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"
                     ),
                     ordered = TRUE)
  ) %>%
  # 物理排序数据框（确保因子顺序与数据行顺序一致）
  arrange(celltype, Dataset)

ggplot(plot_data_clean, aes(x = celltype, y = Value, fill = Dataset)) +  
  geom_col(position = "dodge", width = 0.7) +  # 使用dodge位置并调整柱宽  
  geom_text(aes(label = round(Value, 3)),   
            position = position_dodge(width = 0.7),  
            vjust = -0.5, size = 3) +  
  scale_fill_manual(values = c("#1f78b4", "#33a02c","pink","#f57c6e"),  
                    labels = c("ieu_b_7", 
                               "LARGE_PD_PAPER_TRACTOR_EUR_Meta_analysis_random_LARGE_PD_P1_P2_hg38",
                               "GP2_EUR_ONLY",
                               "Global_Parkinson_Genetics_Program_2023_GWAS_trans_ancestry"
                    )) +  
  labs(title = "sc_ATAC Cell Type Comparison Across GWAS Datasets",  
       x = "Cell Type",  
       y = "-log10(pval)") +  
  theme_minimal() +  
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    legend.title = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  
  ) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 调整y轴范围  
#################################################################



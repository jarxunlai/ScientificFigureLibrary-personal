# =============================================================================
# 堆积柱加折线热图注释
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# =============================================================================

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
set.seed(20260912)
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})

library(ggplot2)
# library(tidyverse) 已在文件头拆包
library(ComplexHeatmap)
library(RColorBrewer)
library(patchwork)

############## 热图数据构建 + 图形绘制 ##################
# 构建模拟数据:
# 上方热图注释：
x = factor(c(paste0("P", 1:17),
             paste0("M", 1:10),
             paste0("C", 1:3)),
             levels = c(paste0("P", 1:17),
                        paste0("M", 1:10),
                        paste0("C", 1:3)))

heatmap_data <- data.frame(
  Data_source = sample(c("Synthetic A", "Synthetic B"), 30, replace = T),
  Sample_Origins = rep(c("Primary", "Distant Metastasis", "Chemotherapy"), c(15,12,3)),
  Smoking = sample(c("Current Smoker", "Former Smoker", "Never Smoker"), 30, replace = T),
  EGFR = sample(c("EGFR Mutation", "WT"), 30, replace = T),
  Stages = sample(paste0("Stage", 1:4), 30, replace = T)
)

rownames(heatmap_data) <- x

cols <- list(Data_source = c("#ff8969", "#e9cd50"),
             Sample_Origins = c("#0e9cc4", "#f3734e", "#c05a9c"),
             Smoking = c("#bed4ad", "#bdd3ac", "#96c18c"),
             EGFR = c("#f7d4b5","#fcf5dd"),
             Stages = c("#feeff4", "#f5cedc", "#f7afcc", "#ec75a7")
)

cols_vec <- unlist(cols)

names(cols_vec) <- c("Synthetic A", "Synthetic B",
                     "Primary", "Distant Metastasis", "Chemotherapy",
                     "Current Smoker", "Former Smoker", "Never Smoker",
                     "EGFR Mutation", "WT",
                     paste0("Stage", 1:4))

# 统一离散横轴顺序，避免热图列和下方样本错位。
heatmap_long <- heatmap_data %>% rownames_to_column("sample") %>% pivot_longer(-sample,names_to="annotation",values_to="value")
heatmap_long$sample <- factor(heatmap_long$sample,levels=levels(x))
p_top <- ggplot(heatmap_long,aes(sample,annotation,fill=value))+geom_tile(colour="white")+scale_fill_manual(values=cols_vec)+theme_minimal()+theme(axis.text.x=element_blank(),axis.title=element_blank(),panel.grid=element_blank(),legend.position="none")

############## 折线图数据构建 + 图形绘制 ##################
# 折线图数据：
line_data <- data.frame(x = factor(c(paste0("P", 1:17),
                                     paste0("M", 1:10),
                                     paste0("C", 1:3)),
                              levels = c(paste0("P", 1:17),
                                         paste0("M", 1:10),
                                         paste0("C", 1:3))),
                        value = c(sort(sample(1:10, 18, replace = T)),
                                  sort(sample(4:12, 9, replace = T)),
                                  sort(sample(6:10, 3, replace = T))
                                  ))

# 折线图：
p2 <- ggplot(line_data, aes(x, value))+
  geom_point(size = 2)+
  geom_line(size = 1, group = 1)+
  scale_y_continuous(breaks = seq(0, 12, 2))+
  ylab("log2(Malignant Cell Number)")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))


############## 堆积柱状图数据构建 + 图形绘制 ##################
# 堆积柱状图数据：
bar_data <- data.frame()

for (i in 1:30) {
  cell_category <- sample(3:6, 1)
  bar_data_tmp <- data.frame(
    x = rep(x[i], cell_category),
    values = sample(20:100, cell_category),
    group = paste0("group", 1:cell_category))
  bar_data_tmp$proportion <- bar_data_tmp$values/sum(bar_data_tmp$values)
  bar_data <- rbind(bar_data, bar_data_tmp)
}

# 堆积柱状图:
p3 <- ggplot(bar_data)+
  geom_bar(aes(x, proportion*100, fill = group),
           color = "white", width = 1, size = 1,
           stat = "identity", position = "stack")+
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  ylab("Cell Proportion(%)")+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))

pdf(file.path(out_dir, "plot.pdf"), height = 6, width = 7)
p_top <- p_top + scale_x_discrete(limits=levels(x),expand=expansion(add=.5))
p2 <- p2 + scale_x_discrete(limits=levels(x),expand=expansion(add=.5))
p3 <- p3 + scale_x_discrete(limits=levels(x),expand=expansion(add=.5)) + theme(axis.text.x=element_text(angle=45,hjust=1))
final_plot <- p_top/p2/p3 + plot_layout(heights=c(1,1.2,2))
print(final_plot)
dev.off()

# AI拼图+调整

ggplot2::ggsave(file.path(root,"preview.png"), plot=final_plot, width=10, height=8, dpi=200, bg="white")

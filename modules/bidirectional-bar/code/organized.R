# =============================================================================
# 双向柱状图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」010双向柱状图
# 本地复现：drafts/ks-shengxin-huitu-repro/010-bidirectional-bar
# Pixi：项目根 default 环境；library(tidyverse) 已拆成 ggplot2/dplyr/tidyr 等。
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

# 构造数据：
data <- data.frame("value" = sample(-5:20, 12))

data$pathway <- paste0("pathway",1:12)

# 初步绘图：
ggplot(data)+
  geom_col(aes(reorder(pathway, value), value))+
  theme_classic()+
  ylim(-20,20)+
  coord_flip()

# 添加颜色，调整主题：
# 设置颜色变量：
color <- rep("#ae4531", 12)
color[which(data$value < 0)] <- "#2f73bb"

data$color <- color


# 最终绘图代码：
ggplot(data)+
  geom_col(aes(reorder(pathway, value), value, fill = color))+
  scale_fill_manual(values = c("#2f73bb","#ae4531"))+
  # 加一条竖线：
  geom_segment(aes(y = 0, yend = 0,x = 0, xend = 12.8))+
  theme_classic()+
  ylim(-20,20)+
  coord_flip()+
  # 调整主题：
  theme(
    # 去除图例：
    legend.position = "none",
    # 标题居中：
    plot.title = element_text(hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
  )+
  ylab("Normalized Enrichment Score")+
  # 添加label：
  geom_text(data=data[which(data$value > 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = 1.1, size = 4)+
  geom_text(data=data[which(data$value < 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = -0.1, size = 4)+
  ggtitle("HER2-enriched Subtype \n FDR < 0.0001")+
  scale_x_discrete(expand=expansion(add=c(0,1.5)))+
  # 添加箭头注释：
  geom_segment(aes(y = -1, yend = -18,x = 13, xend = 13),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  geom_segment(aes(y = 1, yend = 18,x = 13, xend = 13),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  annotate("text", x = 13, y = -20, label = "CP")+
  annotate("text", x = 13, y = 20, label = "MFP")

ggsave(file.path(out_dir, "barplot.pdf"), height = 7, width = 7)

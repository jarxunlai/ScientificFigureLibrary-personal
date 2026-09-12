# =============================================================================
# 渐变柱状图
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

# 加载R包：
library(ggplot2)

# 构建示例数据：
data <- data.frame(
  Gene = paste0("Cancer", 1:21),
  Log2FC = c(runif(15, -0.5,0), runif(6, 0, 0.5)),
  "Log10qvalue" = c(runif(6, 1, 3),runif(15, 0, 1)))

# 查看数据：
head(data)
#      Gene      Log2FC Log10qvalue
# 1 Cancer1 -0.18511148    2.634901
# 2 Cancer2 -0.13808342    1.939038
# 3 Cancer3 -0.30934002    2.632905
# 4 Cancer4 -0.04036049    1.911711
# 5 Cancer5 -0.04722910    1.024237
# 6 Cancer6 -0.18350572    1.759855

ggplot(data)+
  # 竖线：
  geom_hline(yintercept = c(-log2(1.2), -log2(1.5)), 
             color = "grey", linetype = "dotted")+
  # 柱状图
  geom_col(aes(x = reorder(Gene, Log2FC), y = Log2FC, fill = Log10qvalue),
           color = "grey")+
  # 黑线：
  geom_hline(yintercept = 0, color = "black")+
  # 渐变颜色填充：
  scale_fill_gradientn(name="-Log10_\nq-value", # 修改图例标题
                       colours = c("#f6fafd", "#c8dfef", "#6fa6d1", "#2c49a2"),
                       breaks = 0:3,
                       labels = paste0(0:3, ".0"))+
  
  # 坐标轴标题：
  xlab("")+
  ylab("Log2FC(9p21-Loss vs .9p21-WT)")+
  # 增大边缘柱形与边框的距离
  scale_x_discrete(expand = c(0.05, 0.05))+
  # 注释：
  annotate(geom = "text", y = -0.5, x = 21, label = "FC<-1.2")+
  # 标题：
  ggtitle("TCR shannon entroy")+
  # 坐标轴翻转：
  coord_flip()+
  # 主题：
  theme_light()+
  theme(panel.grid = element_blank(),  # 去掉网格线
        plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中、字体
        axis.ticks.y = element_blank(), # 去掉y轴刻度线
        axis.title.x = element_text(size = 10),  # x轴标题大小
        axis.text.y = element_text(size = 10),  # y轴刻度大小
        panel.border = element_rect(color = "black"),
        legend.position=c(0.9,0.1), legend.justification=c(1,0))  # 图例位置

# 保存
ggsave(file.path(out_dir, "barplot.pdf"), height = 5, width = 4)

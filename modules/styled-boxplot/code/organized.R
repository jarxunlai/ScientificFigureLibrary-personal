# =============================================================================
# 美化箱线图
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

library(ggplot2)
# library(tidyverse) 已在文件头拆包
library(ggsignif)

# 数据构建：
# 基础数据矩阵：
data <- data.frame(
  r1 = c(runif(10, min=10, max=40), runif(3, min=0, max=10), runif(3, min=50, max=80)),
  r2 = c(runif(10, min=5, max=30), runif(3, min=0, max=5), runif(3, min=30, max=40)),
  r3 = c(runif(10, min=6, max=30), runif(3, min=0, max=6), runif(3, min=30, max=50)),
  r4 = c(runif(10, min=10, max=40), runif(3, min=0, max=10), runif(3, min=50, max=80)),
  r5 = c(runif(10, min=5, max=30), runif(3, min=0, max=5), runif(3, min=30, max=40)),
  r6 = c(runif(10, min=6, max=30), runif(3, min=0, max=6), runif(3, min=30, max=50)),
  r7 = c(runif(10, min=10, max=100), runif(3, min=0, max=10), runif(3, min=100, max=250)),
  r8 = c(runif(10, min=10, max=40), runif(3, min=0, max=10), runif(3, min=50, max=80)),
  r9 = c(runif(10, min=10, max=40), runif(3, min=0, max=10), runif(3, min=50, max=80)),
  r10 = c(runif(10, min=20, max=150), runif(3, min=0, max=20), runif(3, min=150, max=200)),
  r11 = c(runif(10, min=10, max=50), runif(3, min=0, max=10), runif(3, min=50, max=80)),
  r12 = c(runif(10, min=10, max=40), runif(3, min=0, max=10), runif(3, min=50, max=80))
)

# 宽数据转长数据：
data_long <- data %>%
  pivot_longer(cols = everything(),
               names_to = "group",
               values_to = "values")

data_long$group2 <- factor(rep(rep(c("9p21-WT", "9p21-LOH","9p21-Loss"), 4),16),
                          levels = c("9p21-WT", "9p21-LOH","9p21-Loss"))
data_long$cancer <- rep(rep(c("HNSC_HPV", "PAAD_ALL","SKCM_ALL", "SKCM_BRAF"), each = 3),16)

data_long$group <- factor(data_long$group, levels = paste0("r",1:12))



############# 绘图 -------------
ggplot(data_long)+
  geom_boxplot(aes(group, values, fill = group2, color = group2),
               width = 0.6,
               outlier.shape = NA)+
  geom_line(data = tmp_data, aes(x_value, med_value, group = group),
            color = "white")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(limits = c(0, 240), expand = c(0,0))+
  # 颜色：
  scale_fill_manual(values = c("#bcbdbf", "#86c4e7", "#3675b6"))+
  scale_color_manual(values = c("#bcbdbf", "#86c4e7", "#3675b6"))+
  annotate("rect", xmin = 0.5, xmax = 12.5, ymin = 220, ymax = 240,
           fill = "#eaeae1", color = "#dbd9cd")+
  annotate("text", x = seq(2, 11, 3), y = 230,
           label = unique(data_long$cancer))+
  geom_vline(xintercept = c(3.5, 6.5, 9.5), color = "grey")+
  # 添加p值：
  geom_signif(aes(group, values),
              comparisons = list(c("r1", "r2"),
                                 c("r4", "r5"),
                                 c("r7", "r9"),
                                 c("r10", "r11")),
              vjust = 2,
              tip_length = rep(0,8),
              y_position = c(100, 110, 180, 200),
              # 修改注释内容
              annotations = c("q = 1.5e-5", "q = 2.2e-4",
                              "q = 1.5e-4", "q = 0.008"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none")


# 中位数线数据：
x_value <- c(0.7, 1.3)
for (i in 1:11) {
  x_value <- append(x_value, x_value[c(1, 2)] + i)
}

med_value <- rep(apply(data, 2, median), each = 2)
tmp_data <- data.frame(
  x_value = x_value,
  med_value = med_value,
  group = rep(colnames(data), rep(2,12))
)

# 绘图：
ggplot(data_long)+
  geom_boxplot(aes(group, values, fill = group2, color = group2),
               # 间距调整：
               width = 0.5,
               outlier.shape = NA
               )+
  # 中位数线：
  geom_line(data = tmp_data, aes(x_value, med_value, group = group), color = "#ffffff")+
  # 顶部灰色方块：
  geom_rect(aes(xmin=0, xmax=13, ymin=220, ymax = 240), fill="#eaeae0")+
  # 灰色竖线：
  geom_vline(xintercept = c(3.5, 6.5, 9.5), color = "#bcbdbf", alpha = 0.8)+
  # x轴标签：
  scale_x_discrete(labels = rep(c("9p21-WT", "9p21-LOH","9p21-Loss"), 4))+
  scale_y_continuous(expand = c(0,0))+
  # 颜色：
  scale_fill_manual(values = c("#bcbdbf", "#86c4e7", "#3675b6"))+
  scale_color_manual(values = c("#bcbdbf", "#86c4e7", "#3675b6"))+
  # 主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  # 标题：
  ggtitle("TCR Richness")+
  # 添加p值：
  geom_signif(aes(group, values),
              comparisons = list(c("r1", "r2"),
                             c("r4", "r5"),
                             c("r7", "r9"),
                             c("r10", "r11")),
              vjust = 2,
              tip_length = rep(0,8),
              y_position = c(100, 110, 180, 200),
              # 修改注释内容
              annotations = c("q = 1.5e-5", "q = 2.2e-4",
                              "q = 1.5e-4", "q = 0.008"))+
  annotate("text", x = seq(2,11, 3), y = 230, label = unique(data_long$cancer))


ggsave(file.path(out_dir, "boxplot.pdf"), height = 4, width = 7)

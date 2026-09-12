# =============================================================================
# 多组哑铃图
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

# library(tidyverse) 已在文件头拆包
library(ggplot2)
library(latex2exp)

###### 构造数据：需要p值，log2FC值， 以及分组信息 -----------
# 构造FC值，对应x轴位置：
FC_HFSW <- c(runif(2, -5, -2.5), runif(8, -2.5, 0), runif(10, 0, 5))
FC_SSa <- FC_HFSW + runif(20, -1, 1)
FC_SSd <- FC_HFSW + runif(20, -1, 1)

# 构造P值，对应散点大小：
p_HFSW <- runif(20, 1, 10)*10^(-sample(2:10, 20, replace = T))
p_SSa <- runif(20, 1, 10)*10^(-sample(2:10, 20, replace = T))
p_SSd <- runif(20, 1, 10)*10^(-sample(2:10, 20, replace = T))

# 数据：DAG信息为虚拟构造，无实际意义
data <- data.frame(y = rep(paste0("DAG ", 21:40, ":", 1:20), 3),
                   log2FC = c(FC_HFSW, FC_SSa, FC_SSd),
                   p = c(p_HFSW, p_SSa, p_SSd),
                   group = rep(c("HFSW", "SSa", "SSd"), each = 20))

# 为了绘制横线，需要找出每组中log2FC的最小值和最大值：
data <- data %>%
  group_by(y) %>%
  mutate(x_min = min(log2FC),
         x_max = max(log2FC))

# 查看数据：
head(data)

# # A tibble: 6 × 6
# # Groups:   y [6]
#   y       log2FC         p group  x_min  x_max
#   <chr>      <dbl>         <dbl> <chr>  <dbl>  <dbl>
# 1 DAG 21:1 -4.23   0.00000000724 HFSW  -5.22  -4.00
# 2 DAG 22:2 -3.82   0.000000139   HFSW  -3.82  -2.95
# 3 DAG 23:3 -0.0994 0.0425        HFSW  -0.549  0.615
# 4 DAG 24:4 -1.37   0.00000000403 HFSW  -1.37  -0.392
# 5 DAG 25:5 -1.20   0.0776        HFSW  -2.07  -0.935
# 6 DAG 26:6 -0.770  0.0000224     HFSW  -1.73  -0.770

# 保存数据：

########## 绘图 -----------
ggplot(data)+
  # 散点图：
  geom_point(aes(log2FC, y, fill = group, size = -log10(p)),
             shape = 21) +
  # 横线：
  geom_linerange(aes(xmin = x_min, xmax = x_max,
                     y = y), linewidth = 0.3)+
  # 竖虚线：
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey", linewidth = 0.5)+
  # 散点填充色：
  scale_fill_manual(name = "",
                    values = c("#f5b4a9", "#b2c6e6", "#7abc7c"))+
  # 大小设置，及图例设置：
  scale_size_continuous(name = expression(-log[10](P)),
                        range = c(1,3))+
  # 坐标轴标签：
  ylab("")+
  xlab(TeX("$Log_{2}FC$", bold = T))+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(face = "bold"))

# 保存：
ggsave(file.path(out_dir, "plot.pdf"), height = 6, width = 5)

ggplot2::ggsave(file.path(root, "preview.png"), plot = ggplot2::last_plot(), width = 5, height = 6, dpi = 200, bg = "white")

# =============================================================================
# 散点饼图
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

library(scatterpie)
# library(tidyverse) 已在文件头拆包
library(readxl)

# 独立模拟的组成数据，所有成分每行合计 100%。
data <- read.csv(file.path(root,"data","input.csv"), check.names=FALSE)

# 转换x轴信息为因子，最终转为数值：
data$Time <- factor(data$Time, levels = c("d1", "d14", "d30", "d60", "d90", ">d180"))
data <- data[order(data$Time), ]
data$x <- 1:nrow(data)-0.5
data$region <- factor(1:nrow(data))

# 去掉前两列：
data_new <- data



# 设置颜色：
col_values <- c("#325939", "#6eaa62", "#b4d79f", "#f2d533",
                "#8daed1", "#6752a1", "#9d9ac5", "#e36944")
# 绘图：
ggplot()+
  # 散点饼图：
  geom_scatterpie(aes(x = x, y = Proportion, group = region),
                  color = NA, pie_scale = 1,
                  cols = colnames(data_new)[3:10],
                  data = data_new) +
  # 虚线：
  geom_vline(xintercept = cumsum(as.numeric(table(data_new$Time)))[-length(table(data_new$Time))],
             linetype = "dashed", color = "grey")+
  # x轴刻度调节：
  scale_x_continuous(breaks = c(2, 6, 10.5, 15, 19.5, 23),
                     labels = c("d1\n(n = 4)", "d14\n(n = 4)", "d30\n(n = 5)",
                                "d60\n(n = 4)", "d90\n(n = 5)", ">d180\n(n = 2)"),
                     expand = c(0,0))+
  # 主题和标签：
  ylab("Simulated proportion (%)")+
  scale_fill_manual(name = "Cell type", values = col_values)+
  coord_equal()+
  theme_classic()

# 保存
ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 7)

ggplot2::ggsave(file.path(root,"preview.png"), plot=ggplot2::last_plot(), width=7, height=5, dpi=200, bg="white")

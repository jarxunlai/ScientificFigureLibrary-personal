# =============================================================================
# 冲击图
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
library(ggforce)
# library(tidyverse) 已在文件头拆包

# 构造模拟数据：
data <- data.frame(Normal = sample(1:100, 8),
                   A = sample(1:100, 8),
                   B = sample(1:100, 8),
                   C = sample(1:100, 8),
                   group = factor(c("Epithelial", "Endothelial", "Myofibroblast",
                             "Fibroblast", "Myeloid", "T_NK_Cell", "B_cell",
                             "Plasma_cell"),
                             levels = c("Epithelial", "Endothelial", "Myofibroblast",
                                        "Fibroblast", "Myeloid", "T_NK_Cell", "B_cell",
                                         "Plasma_cell")))

# 宽数据转长数据：
data_long <- pivot_longer(data, cols = !group,
                          names_to = "x", values_to = "value") %>%
  group_by(x)
data_long$x <- factor(data_long$x, levels = c("Normal", "A", "B", "C"))

# 构造过度条带数据：
data_y_tmp <- data[8:1,]

for (j in 1:4) {
  for(i in 2:nrow(data_y_tmp)) {
    data_y_tmp[i, j] <- sum(data[8:1,][1:i,j])
  }
}
data_y_tmp[,1:4] <- apply(data_y_tmp[,1:4], 2, function(x) x/x[8])

y <- c()

for (j in 1:3) {
  for (i in 1:8) {
    if (i == 1) {
      y[((i-1)*4+(j-1)*32+1):(i*4+(j-1)*32)] <- c(0, 0, data_y_tmp[i,j+1], data_y_tmp[i,j])
    } else {
      y[((i-1)*4+(j-1)*32+1):(i*4+(j-1)*32)] <- c(data_y_tmp[i-1,j], data_y_tmp[i-1,j+1], data_y_tmp[i,j+1], data_y_tmp[i,j])
    }
  }
}

data_y <- data.frame(x = rep(rep(c(1.25, 1.75, 1.75, 1.25), 8), 3) + rep(0:2, each = 32),
                     group1 = rep(1:24, each = 4),
                     group2 = rep(rep(factor(rev(c("Epithelial", "Endothelial", "Myofibroblast",
                                               "Fibroblast", "Myeloid", "T_NK_Cell", "B_cell",
                                               "Plasma_cell")),
                                             levels = c("Epithelial", "Endothelial", "Myofibroblast",
                                                        "Fibroblast", "Myeloid", "T_NK_Cell", "B_cell",
                                                        "Plasma_cell")), each = 4), 3),
                     y = y)


# 绘图：
ggplot()+
  geom_col(data = data_long, aes(x, value, fill = group),
           position = "fill", width = 0.5)+

  geom_diagonal_wide(data = data_y,
                     aes(x, y, group = group1, fill = group2),
                     alpha = 0.4)+
  scale_x_discrete(label = c("Normal", "0+", "1+", "2+"))+
  scale_fill_manual(name = "Annol",
                    values = c("#e9cdb4", "#016a99", "#d49d4d","#abddde",
                               "#fd0103", "#029f88", "#f4ac05", "#f58502"))+
  ylab("Propotion")+
  xlab("HER2")+
  theme_classic()

ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 7)

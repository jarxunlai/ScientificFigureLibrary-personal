# =============================================================================
# Upset图
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

library(UpSetR)
# 载入R包：
require(ggplot2)
require(plyr)
require(gridExtra)
require(grid)

# 载入示例电影类型表（从 UpSetR extdata 拷入 data/）：
movies <- read.csv(file.path(root, "data", "movies.csv"), header = TRUE, sep = ";")

########## 基础绘图 ----------------
pdf(file.path(out_dir, "plot1.pdf"), height = 5, width = 7)
upset(movies,
      # 数据集数量：
      nsets = 7,
      # 柱形的最大数目：
      nintersects = 30,
      # 柱形图与矩阵图大小比例：
      mb.ratio = c(0.5, 0.5),
      # 柱形的排序方式：freq -- 柱形高度排序； degree -- 数据集的数量排序；
      order.by = c("freq", "degree"), decreasing = c(TRUE, FALSE))
dev.off()

########## 美化 -------------------
pdf(file.path(out_dir, "plot2.pdf"), height = 5, width = 7)
upset(movies,
      # 数据集数量：
      nsets = 5,
      # 柱形图与矩阵图大小比例：
      mb.ratio = c(0.6, 0.4),
      # 修改左侧柱形颜色：
      sets.bar.color = c("#006bba", "#d4550e", "#00a6dc", "#67a330", "#f88a53"),
      # 柱形排序：
      order.by = "freq", decreasing = T,
      # 修改条形的颜色和矩阵散点的颜色：
      queries = list(list(query = intersects, params = list("Drama"),
                          color= "#006bba", active = T)),
      # 矩形散点的阴影颜色：
      shade.color = NA
      )
dev.off()

######### 批量修改 ---------------
# 先创建颜色，有多少个柱子就设置多少个颜色：
library(RColorBrewer)

colors <- colorRampPalette(brewer.pal(9, "Set1"))(27)

query_list <- list()

tmp <- unique(movies[, c("Drama", "Comedy", "Action", "Thriller", "Romance")])
tmp <- tmp[rowSums(tmp) != 0, ]

for (i in 1:27) {
  query_list[[i]] <- list(query = intersects,
                          params = list(colnames(tmp)[which(tmp[i,] == 1)]),
                          color= colors[i], active = T)
}

pdf(file.path(out_dir, "plot.pdf"), height = 3, width = 6)
upset(movies,
      # 数据集数量：
      nsets = 5,
      # 柱形图与矩阵图大小比例：
      mb.ratio = c(0.65, 0.35),
      # 修改左侧柱形颜色：
      sets.bar.color = colors[sample(1:27, 5)],
      # 柱形排序：
      order.by = "freq", decreasing = T,
      # 修改条形的颜色和矩阵散点的颜色：
      queries = query_list,
      # 矩形散点的阴影颜色：
      shade.color = NA
)
dev.off()

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
if (exists("movies")) {
  png(file.path(out_dir, "upset1.png"), width = 1400, height = 900, res = 150)
  print(UpSetR::upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE, FALSE)))
  dev.off()
}

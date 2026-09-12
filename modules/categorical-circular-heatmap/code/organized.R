# =============================================================================
# 分类变量环形热图
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

# 载入R包：
library(readxl)
library(circlize)
# library(tidyverse) 已在文件头拆包

data <- read.csv(file.path(root, "data", "input.csv"), row.names = 1, check.names = FALSE)

# 设置颜色模式：
col_fun = list(col_1 = c("Y" = "#ee5567", "N" = "#fbeaeb"),
               col_2 = c("Y" = "#fc6f53", "N" = "#faeee7"),
               col_3 = c("Y" = "#fdcb55", "N" = "#fdf7e8"),
               col_4 = c("Y" = "#9ed365", "N" = "#f2f9ea"),
               col_5 = c("Y" = "#4bdde7", "N" = "#e6f5fa"),
               col_6 = c("Y" = "#5e99e6", "N" = "#ebf3f9"),
               col_7 = c("Y" = "#7d66b7", "N" = "#eeedf4")
)


# 绘图：
{
  grDevices::png(file.path(root, "preview.png"), width = 1800, height = 1800, res = 150, type = "cairo")
  circos.par(gap.degree = 10, start.degree = 85, track.margin = c(0.001, 0.001))

  data_tmp <- as.matrix(data[,8])
  rownames(data_tmp) <- rownames(data)
  circos.heatmap.initialize(data_tmp, cluster = F)

  circos.track(
    ylim = c(0, 6),
    bg.border = NA,
    panel.fun = function(x, y){
      circos.barplot(data_tmp,
                     border = NA,
                     1:nrow(data_tmp) - 0.5,
                     col = "#b7a085")
      })

for (i in 1:7) {
  data_tmp <- as.matrix(data[,i])
  if (i == 1) {
    rownames(data_tmp) <- rownames(data)
  }
  colnames(data_tmp) <- colnames(data)[i]

  circos.heatmap(data_tmp,
                 col = col_fun[[i]],
                 rownames.side = "outside",
                 cluster = FALSE,
                 cell.border = "white",
                 track.height = 0.05)
}

circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                CELL_META$cell.xlim[2] + convert_x(6, "mm"), 7,
                col = "#e5e8ed", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 0.5,
                "7", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 1.5,
                "6", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 2.5,
                "5", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 3.5,
                "4", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 4.5,
                "3", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 5.5,
                "2", cex = 0.5, facing = "inside")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), 6.5,
                "1", cex = 0.5, facing = "inside")

  }
}, bg.border = NA)
circos.clear()
dev.off()
}

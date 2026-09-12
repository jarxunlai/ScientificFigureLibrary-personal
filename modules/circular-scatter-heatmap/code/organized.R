# =============================================================================
# 环形散点加热图
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

library(readxl)
# library(tidyverse) 已在文件头拆包

data <- read.csv(file.path(root, "data", "input.csv"), check.names = FALSE)

# 绘图:
library(circlize)

################### 外圈热图 ========================
# 设置颜色模式:
col_fun = list(col_1 = colorRamp2(c(0, 0.5, 1),
                                  c("#ac90f4", "#8df4e6", "#e99651")),
               col_2 = colorRamp2(c(0, 1),
                                  c("#037b3d", "#55cbd7")),
               col_3 = colorRamp2(c(0, 0.5, 1),
                                  c("#15489c", "#fff403", "#e93832")),
               col_4 = colorRamp2(c(0, 1),
                                  c("#3c9d9d", "#cacaff")),
               col_5 = colorRamp2(c(0, 0.5, 1),
                                  c("#fe34cb", "#c0b892", "#6b396b"))
)


if (T) {
  grDevices::png(file.path(root, "preview.png"), width = 1800, height = 1800, res = 150, type = "cairo")
  ################### 外圈热图 ========================
  circos.par("track.height" = 0.1,   # 每行的高度
             "start.degree" = 0,  # 环形开始的角度
             "gap.degree" = c(5, 5, 5, 5, 30), # sector之间的gap角度
             "track.margin" = c(0.01, 0.01))  # 每行之间的距离

  data_1 <- data[,-4]

  # 循环绘制热图：
  for (i in 1:6) {
    data_tmp <- as.matrix(as.numeric(data_1[, i+2]))
    if (i == 1) {
      rownames(data_tmp) <- data$ID
    }
    colnames(data_tmp) <- paste0(colnames(data_tmp), "-p")

    if (i < 6) {
      circos.heatmap(data_tmp,
                     # 设置sector：
                     split = data$Level,
                     col = col_fun[[i]],
                     rownames.side = "outside",
                     rownames.cex = 0.3,  # 行名字体大小
                     cluster = F,
                     cell.border = NA,  # 单元格边框
                     track.height = 0.05)
    } else {
      ################### 内圈散点图 ========================
      data_tmp <- as.matrix(as.numeric(data[, 4]))
      colnames(data_tmp) <- "IVW-OR"

      circos.track(ylim = range(data_tmp[,1]),
                   panel.fun = function(x, y) {
                     # sector的名称：
                     circos.text(CELL_META$xcenter,
                                 CELL_META$cell.ylim[2] + mm_y(1.5),
                                 CELL_META$sector.index,
                                 cex = 0.5
                     )
                     # y轴刻度：
                     circos.yaxis(labels.cex = 0.3,
                                  at = seq(0.6, 1.8, 0.4),
                                  sector.index = "phylum")
                     y = data_tmp[,1][CELL_META$subset]
                     y = y[CELL_META$row_order]
                     # 虚线：
                     circos.lines(CELL_META$cell.xlim, c(1, 1), lty = "dashed", col = "black")
                     # 散点：
                     circos.points(seq_along(y) - 0.5, y,
                                   pch = 16, cex = 0.5,
                                   col = "#63caff")
                     # 行名：
                     if(CELL_META$sector.numeric.index == 5) { # the last sector
                       cn = rev(c("IVW-p", "MR Egger-p", "WM-p", "MLE-p","MR RAPS-p", "IVW-OR"))
                       n = length(cn)
                       circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(8, "mm"),
                                   c(1, seq(3.5, 9.5, 1.5)), cn,
                                   cex = 0.5, adj = 0.5, facing = "inside")
                     }
                   },
                   track.margin = c(0.02, 0.02),
                   cell.padding = c(0.02, 0, 0.02, 0))
    }
  }
  circos.clear()
  dev.off()
}

# =============================================================================
# 复杂环形热图
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

data_top200 <- read.csv(file.path(root, "data", "input.csv"), row.names = 1, check.names = FALSE)

# 绘图：
library(circlize)

col_7 <- c("#a6cee3", "#1f78b4", "#b2df8a",
           "#33a02c", "#fb9a99", "#e31a1c",
           "#fdbf6f", "#ff7f00", "#cab2d6",
           "#6a3d9a")
names(col_7) <- rev(unique(data_top200$SUPER_PATHWAY))


col_fun = list(col_1 = colorRamp2(c(-11, 0, 11),
                                  c("#0775b4", "white", "#a54b05")),
               col_2 = colorRamp2(c(-4, 0, 6),
                                  c("#840822", "white", "#181818")),
               col_3 = colorRamp2(c(0, 0.15, 0.3),
                                  c("white", "#41ab5d", "#00431a")),
               col_4 = colorRamp2(c(0, 0.15, 0.3),
                                  c("white", "#f5563c", "#68000d")),
               col_5 = colorRamp2(c(0, 0.15, 0.3),
                                  c("white", "#9a96c6", "#3f007d")),
               col_6 = colorRamp2(c(0, 0.15, 0.3),
                                  c("white", "#4f9acb", "#08306b")),
               col_7 = col_7
               )
if (T) {
  grDevices::png(file.path(root, "preview.png"), width = 1800, height = 1800, res = 150, type = "cairo")
  circos.par(gap.degree = 70, start.degree = 20, track.margin = c(0.001, 0.001))

  for (i in 1:7) {
    data_tmp <- as.matrix(data_top200[,i])
    if (i == 1) {
      rownames(data_tmp) <- rownames(data_top200)
    }
    colnames(data_tmp) <- colnames(data_top200)[i]

    if (i < 7) {
      circos.heatmap(data_tmp,
                     col = col_fun[[i]],
                     rownames.side = "outside",
                     cluster = F,
                     cell.border = "white",
                     track.height = 0.05)
    } else {
      circos.heatmap(data_tmp,
                     col = col_fun[[i]],
                     rownames.side = "outside",
                     cluster = F,
                     track.height = 0.03)
    }
  }

  # 图例：
  library(ComplexHeatmap)

  # Non ACS versus ACS-log10(P)
  lgd1 <- Legend(title = "", border = "black", grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(-11, 0, 11), title_position = "topcenter",
                 col_fun = col_fun[[1]], direction = "horizontal")

  lgd2 <- Legend(title = "", border = "black", grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(-4, 0, 6), title_position = "topcenter",
                 col_fun = col_fun[[2]], direction = "horizontal")

  lgd3 <- Legend(title = "", border = "black", grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(0, 0.3), title_position = "topcenter",
                 col_fun = col_fun[[3]], direction = "horizontal")

  lgd4 <- Legend(title = "", border = "black",
                 grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(0, 0.3), title_position = "topcenter",
                 col_fun = col_fun[[4]], direction = "horizontal")

  lgd5 <- Legend(title = "", border = "black", grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(0, 0.3), title_position = "topcenter",
                 col_fun = col_fun[[5]], direction = "horizontal")

  lgd6 <- Legend(title = "", border = "black", grid_height = unit(3, "mm"),
                 legend_width = unit(20, "mm"),
                 at = c(0, 0.3), title_position = "topcenter",
                 col_fun = col_fun[[6]], direction = "horizontal")

  pd <- packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, lgd6, row_gap = unit(1, "mm"))
  draw(pd, x = unit(0.55, "npc"), y = unit(0.7, "npc"))

  lgd = Legend(labels = names(col_7),
               legend_gp = gpar(fill = col_7), title = "Super pathways",
               ncol = 1, row_gap = unit(1, "mm"))

  draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))

  dev.off()
  circos.clear()
}

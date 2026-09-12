# =============================================================================
# 三维嵌入坐标散点图（模拟）
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

library(scatterplot3d)
# 三条模拟分支坐标；未运行 PCA 或扩散映射算法。
t <- rep(seq(0,1,length.out=100),3)
group <- rep(1:3,each=100)
plot.data <- data.frame(DC_1=t*cos((group-1)*2*pi/3)+rnorm(300,0,.04), DC_2=t*sin((group-1)*2*pi/3)+rnorm(300,0,.04), DC_3=.3*t+rnorm(300,0,.04), color=c("#FD6AB0","#67C5E8","#87C875")[group])

# scatterplot3d包：
png(file.path(root,"preview.png"),width=1400,height=1400,res=200,type="cairo")
scatterplot3d(x = plot.data$DC_1,
              y = plot.data$DC_3,
              z = plot.data$DC_2,
              color = plot.data$color,
              pch = 16, cex.symbols = 1,
              scale.y = 0.7, angle = 120,
              xlab = "DC_1", ylab = "DC_3", zlab = "DC_2",
              col.axis = "#444444", col.grid = "#CCCCCC")
dev.off()

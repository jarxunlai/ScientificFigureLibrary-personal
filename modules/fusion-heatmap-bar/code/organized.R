# =============================================================================
# 基因融合热图堆积柱
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


library(ComplexHeatmap)
# library(tidyverse) 已在文件头拆包
library(ggplot2)
library(circlize)
library(patchwork)


# 仅发布按基因和肿瘤类型聚合的计数矩阵，不包含个体记录。
data_mat <- read.csv(file.path(root,"data","input.csv"),row.names=1,check.names=FALSE)

data_mat <- data_mat[,order(colSums(data_mat), decreasing = T)]
data_mat <- data_mat[order(rowSums(data_mat), decreasing = T),]
data_mat2 <- as.matrix(data_mat)
data_mat2[which(data_mat2 == 0)] <- NA

# 基础热图：
Heatmap(data_mat2,
        na_col = "white",
        # 去掉行列聚类：
        cluster_rows = F,
        cluster_columns = F)


# 修改颜色+添加文字和描边
col_fun = colorRamp2(c(0, 5, 30, max(data_mat)), c("#b4d9e5", "#91a1cf", "#716bbf","#5239a3"))

p1 <- Heatmap(data_mat2,
        col = col_fun,
        na_col = "white",
        # 去掉行列聚类：
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left", show_row_names = FALSE,
        # 图例
        heatmap_legend_param = list(
          title = "Count",
          title_position = "leftcenter",
          legend_direction = "horizontal"
        ),
        # 行名和列名：
        row_names_gp = gpar(fontsize = 10, font = 3),
        column_names_gp = gpar(fontsize = 10, font = 3),
        # 添加文字注释：
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (!is.na(data_mat2[i,j])) {
            grid.text(sprintf("%1.f", data_mat2[i, j]), x, y,
                      gp = gpar(fontsize = 10, col = "#df9536"))
            grid.rect(x, y, width, height,
                      gp = gpar(col = "grey", fill = NA, lwd = 0.8))
          }
        })

pdf(file.path(out_dir, "Heatmap.pdf"), height = 8, width = 8)
draw(p1, heatmap_legend_side = "bottom")
dev.off()


# 同一个行顺序上的组成堆积柱，与热图逐行对齐。
composition <- as.matrix(data_mat/rowSums(data_mat))
tumor_colors <- setNames(grDevices::hcl.colors(ncol(composition),"Dynamic"),colnames(composition))
left_bar <- rowAnnotation(Gene=anno_text(rownames(data_mat),just="right",location=unit(1,"npc"),gp=gpar(fontsize=10,fontface="italic")),Composition=anno_barplot(composition,gp=gpar(fill=tumor_colors,col=NA),bar_width=.8,width=unit(3,"cm")))
composition_legend <- Legend(title="Tumor type",labels=names(tumor_colors),legend_gp=gpar(fill=tumor_colors),ncol=3)
png(file.path(root,"preview.png"),width=2000,height=1900,res=200,type="cairo")
draw(left_bar+p1,heatmap_legend_side="bottom",annotation_legend_side="bottom",annotation_legend_list=list(composition_legend))
dev.off()

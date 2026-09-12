# =============================================================================
# 热图局部放大
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

# 加载R包：
library(ComplexHeatmap)
library(grid)

# 完全模拟的突变矩阵；静态总览和局部放大，不启动交互服务。
mat <- matrix(sample(c("","HOMDEL","AMP","MUT"),26*172,TRUE,prob=c(.9,.025,.025,.05)),nrow=26,dimnames=list(paste0("Gene_",1:26),paste0("sim_",1:172)))

col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33,
              gp = gpar(fill = col["MUT"], col = NA))
  }
)

column_title = "Simulated alterations: overview"
heatmap_legend_param = list(title = "Alterations", at = c("HOMDEL", "AMP", "MUT"),
                            labels = c("Deep deletion", "Amplification", "Mutation"))
ht = oncoPrint(mat,
               alter_fun = alter_fun, col = col, row_order=seq_len(nrow(mat)), column_order=seq_len(ncol(mat)),
               remove_empty_columns = TRUE, remove_empty_rows = TRUE,
               top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                  foo1 = 1:172,
                                                  bar1 = anno_points(1:172)
               ),
               left_annotation = rowAnnotation(foo2 = 1:26),
               right_annotation = rowAnnotation(bar2 = anno_barplot(1:26)),
               column_title = column_title, heatmap_legend_param = heatmap_legend_param)

# 相同基因顺序的前 20 个模拟样本局部放大。
zoom <- oncoPrint(mat[,1:20],alter_fun=alter_fun,col=col,row_order=seq_len(nrow(mat)),column_order=1:20,show_column_names=FALSE,column_title="Zoom: samples 1-20",heatmap_legend_param=heatmap_legend_param)
png(file.path(root,"preview.png"),width=2400,height=1500,res=180,type="cairo")
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2,widths=unit(c(.65,.35),"npc"))))
pushViewport(viewport(layout.pos.col=1)); draw(ht,newpage=FALSE); popViewport()
pushViewport(viewport(layout.pos.col=2)); draw(zoom,newpage=FALSE,show_heatmap_legend=FALSE); popViewport(2)
dev.off()

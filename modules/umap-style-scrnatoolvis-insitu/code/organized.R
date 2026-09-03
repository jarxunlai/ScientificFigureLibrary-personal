# =============================================================================
# 簇内直接标注 UMAP（scRNAtoolVis 风格）
# organized 版本：中文分节导航 + 线性脚本，不封装函数
# =============================================================================
# 与 code/original.R 的关系：
#   original.R 是微信文章公开的作者脚本，依赖 UMAP.Rdata / Seurat / 专用可视化包。
#   本整理版只保留绘图层：读取合成 UMAP 坐标表，不加载单细胞对象。
# 数据说明：
#   细胞数与图 3/4 图例一致（n=45300）。坐标是按截图簇位生成的合成高斯点云，
#   不是原文 UMAP.Rdata。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- 第一步：定位条目根目录 ------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) %in% c("code", "validation")) {
  dirname(script_dir)
} else {
  script_dir
}

# ---- 第二步：读取合成 UMAP 坐标 --------------------------------------------
plot_df <- read.csv(file.path(root, "data/example_umap.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("UMAP_1", "UMAP_2", "celltype") %in% names(plot_df)))

# ---- 第三步：按截图配色，簇中心放标签 ------------------------------------
# 原文 clusterCornerAxes 依赖 scRNAtoolVis；整理版用 ggplot 近似：点 + 簇内标签 + 角落坐标轴。
cell_levels <- c(
  "Fibroblasts", "T-cells", "Pericyte", "NK-cells", "B-cells",
  "myeloid cells", "Endothelial cells", "smooth muscle cells", "cardiomyocyte cells"
)
pal <- c(
  "Fibroblasts"         = "#1F77B4",
  "T-cells"             = "#2CA02C",
  "Pericyte"            = "#D62728",
  "NK-cells"            = "#9467BD",
  "B-cells"             = "#E377C2",
  "myeloid cells"       = "#BCBD22",
  "Endothelial cells"   = "#17BECF",
  "smooth muscle cells" = "#AEC7E8",
  "cardiomyocyte cells" = "#FFBB78"
)
plot_df$celltype <- factor(plot_df$celltype, levels = cell_levels)
label_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# ---- 第四步：角落坐标轴 UMAP ----------------------------------------------
p1 <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 0.15, alpha = 0.7, stroke = 0) +
  geom_text(
    data = label_df, aes(label = celltype),
    size = 3.4, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = pal, name = "celltype") +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  coord_equal() +
  annotate("segment", x = -12, xend = -7.5, y = -12, yend = -12,
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
  annotate("segment", x = -12, xend = -12, y = -12, yend = -7.5,
           arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
  annotate("text", x = -9.7, y = -12.7, label = "UMAP1", size = 3.2) +
  annotate("text", x = -13.1, y = -9.7, label = "UMAP2", size = 3.2, angle = 90) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.margin = margin(8, 8, 18, 28)
  )

# ---- 最后一步：导出 PNG / PDF ----------------------------------------------
out_dir <- file.path(root, "validation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
png_path <- file.path(out_dir, "synthetic_render.png")
pdf_path <- file.path(out_dir, "synthetic_render.pdf")
preview_path <- file.path(root, "preview.png")
ggsave(filename = png_path, plot = p1, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(filename = pdf_path, plot = p1, width = 8, height = 6, bg = "white")
file.copy(png_path, preview_path, overwrite = TRUE)
message("wrote ", png_path)

# =============================================================================
# 黑边颗粒风 UMAP（SCpubr 风格）
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
  library(shadowtext)
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

# ---- 第三步：shape=21 黑边 + 簇中心描边标签 --------------------------------
# 原文还调用了未定义的 label_df；整理版在这里现算中位坐标。
cell_levels <- c(
  "Fibroblasts", "T-cells", "Pericyte", "NK-cells", "B-cells",
  "myeloid cells", "Endothelial cells", "smooth muscle cells", "cardiomyocyte cells"
)
pal <- c(
  "Fibroblasts"         = "#1F77B4",
  "T-cells"             = "#FF7F0E",
  "Pericyte"            = "#2CA02C",
  "NK-cells"            = "#D62728",
  "B-cells"             = "#9467BD",
  "myeloid cells"       = "#8C564B",
  "Endothelial cells"   = "#E377C2",
  "smooth muscle cells" = "#7F7F7F",
  "cardiomyocyte cells" = "#BCBD22"
)
plot_df$celltype <- factor(plot_df$celltype, levels = cell_levels)
label_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, fill = celltype)) +
  geom_point(shape = 21, color = "black", stroke = 0.12, size = 0.6, alpha = 0.85) +
  shadowtext::geom_shadowtext(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = celltype),
    color = "black", bg.color = "white", size = 5, fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = pal, name = "Cell type") +
  labs(title = "Celltype") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3)))

# ---- 最后一步：导出 PNG / PDF ----------------------------------------------
out_dir <- file.path(root, "validation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
png_path <- file.path(out_dir, "synthetic_render.png")
pdf_path <- file.path(out_dir, "synthetic_render.pdf")
preview_path <- file.path(root, "preview.png")
ggsave(filename = png_path, plot = p1, width = 8, height = 7, dpi = 300, bg = "white")
ggsave(filename = pdf_path, plot = p1, width = 8, height = 7, bg = "white")
file.copy(png_path, preview_path, overwrite = TRUE)
message("wrote ", png_path)

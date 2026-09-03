# =============================================================================
# 半透明椭圆加同色标签 UMAP
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
  library(ggrepel)
  library(colorspace)
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

# ---- 第三步：Dark 3 配色、椭圆与中位标签 ----------------------------------
plot_df$celltype <- factor(plot_df$celltype)
celltype_levels <- levels(plot_df$celltype)
celltype_colors <- setNames(
  qualitative_hcl(length(celltype_levels), palette = "Dark 3"),
  celltype_levels
)
label_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p1 <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = celltype, fill = celltype)) +
  stat_ellipse(geom = "polygon", level = 0.65, alpha = 0.12, color = NA, show.legend = FALSE) +
  geom_point(size = 0.25, alpha = 0.35, stroke = 0, shape = 16) +
  ggrepel::geom_text_repel(
    data = label_df, aes(label = celltype),
    size = 5, fontface = "bold", segment.color = NA,
    max.overlaps = Inf, show.legend = FALSE, seed = 42
  ) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) +
  labs(title = "Celltype", x = "UMAP_1", y = "UMAP_2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 14),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_cartesian(clip = "off")

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

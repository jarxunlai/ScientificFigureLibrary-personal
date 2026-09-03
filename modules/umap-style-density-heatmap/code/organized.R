# =============================================================================
# 暗夜密度热力 UMAP
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
  library(viridis)
  library(ggrepel)
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

# ---- 第三步：magma 密度 + 白色虚线椭圆 ------------------------------------
umap_data <- plot_df %>% rename(cell_type = celltype)
label_centroids <- umap_data %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
x_range <- range(umap_data$UMAP_1)
y_range <- range(umap_data$UMAP_2)
pad <- 0.05 * diff(x_range)

p1 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(geom = "raster", aes(fill = after_stat(density)), contour = FALSE, n = 200) +
  geom_point(color = "#FFFFFF33", size = 0.03, alpha = 0.15) +
  stat_ellipse(
    aes(group = cell_type),
    color = "white", linetype = "dashed", alpha = 0.5, linewidth = 0.55
  ) +
  geom_text_repel(
    data = label_centroids, aes(label = cell_type),
    color = "white", size = 5, fontface = "bold",
    box.padding = unit(0.6, "lines"), max.overlaps = Inf,
    force = 8, segment.color = "white", segment.alpha = 0.5
  ) +
  scale_fill_viridis(option = "magma", direction = 1) +
  coord_cartesian(xlim = x_range + c(-pad, pad), ylim = y_range + c(-pad, pad), expand = TRUE) +
  theme_void(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10, "pt")
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

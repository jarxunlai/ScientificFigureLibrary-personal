# =============================================================================
# 编号加侧边图例 UMAP（SCP 风格）
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

# ---- 第三步：编号标签 + 带细胞数的图例 ------------------------------------
# 原文 SCP::CellDimPlot；整理版用 ggplot 编号与图例复刻信息密度。
cell_levels <- c(
  "Fibroblasts", "T-cells", "Pericyte", "NK-cells", "B-cells",
  "myeloid cells", "Endothelial cells", "smooth muscle cells", "cardiomyocyte cells"
)
pal <- c(
  "Fibroblasts"         = "#AEC7E8",
  "T-cells"             = "#1F77B4",
  "Pericyte"            = "#98DF8A",
  "NK-cells"            = "#2CA02C",
  "B-cells"             = "#FFBB78",
  "myeloid cells"       = "#FF7F0E",
  "Endothelial cells"   = "#E377C2",
  "smooth muscle cells" = "#D62728",
  "cardiomyocyte cells" = "#C5B0D5"
)
plot_df$celltype <- factor(plot_df$celltype, levels = cell_levels)
n_tab <- plot_df %>% count(celltype, name = "n")
legend_lab <- paste0(as.integer(n_tab$celltype), ": ", n_tab$celltype, "(", n_tab$n, ")")
names(legend_lab) <- as.character(n_tab$celltype)
label_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop") %>%
  mutate(idx = as.integer(celltype))

p1 <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 0.18, alpha = 0.8, stroke = 0) +
  geom_text(data = label_df, aes(label = idx), color = "grey15", size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = pal, labels = legend_lab, name = "celltype:") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = paste0("nCells:", nrow(plot_df)), x = "umap_1", y = "umap_2") +
  coord_equal() +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 12),
    legend.text = element_text(size = 9)
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

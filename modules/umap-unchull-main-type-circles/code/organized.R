# =============================================================================
# 单细胞 UMAP 大群虚线非凸包 + 亚群数字标签
# organized 版本：中文分节导航，线性脚本，不新增绘图封装函数
# =============================================================================
# 来源：微信公众号《高分杂志同款单细胞umap加圈：升级版本》
#   生信技能树 / 新叶910，2026-04-27
#   https://mp.weixin.qq.com/s/npmNC0xF8juoxc_xoo6CiA
# 文献：Liu et al. Cell Stem Cell 2025, Figure 1B
#   Spatiotemporal single-cell roadmap of human skin wound healing
#   GEO: GSE241132（58,823 细胞；本仓库没有该矩阵）
# 与 original.R 的关系：
#   original.R 是网页代码原样，依赖 Seurat / 10X / 作者注释表。
#   本整理版只保留绘图层：读取合成 UMAP 坐标，用 ggunchull::stat_unchull
#   和 tidydr::theme_dr。坐标按教程终图簇位生成，不是原文降维。
# 网页 names(cols1) 把 Myeloid 的 22 写成了 2，与 2 号成纤维重复；
# 这里按注释改回 22。Hair follicle 出现在教程终图标签，不在网页 cols2，
# 另加一色。
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggunchull)
library(tidydr)
library(grid)

# ---- 路径 ------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) %in% c("code", "validation")) {
  dirname(script_dir)
} else if (file.exists(file.path(script_dir, "data", "example_umap.csv"))) {
  script_dir
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

# ---- 读取合成坐标 ----------------------------------------------------------
df <- read.csv(
  file.path(root, "data", "example_umap.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(all(c("umap_1", "umap_2", "RNA_snn_res.0.5", "newMainCellTypes") %in% names(df)))
df$RNA_snn_res.0.5 <- factor(df$RNA_snn_res.0.5)

# ---- 亚群色 / 大群色（文献抠色，修正 22） --------------------------------
cols1 <- c(
  "0" = "#c04e2a", "1" = "#f3b169", "3" = "#f0924f", "7" = "#f5c28c",
  "9" = "#cfb8cf", "11" = "#b4ce8e", "14" = "#3d8b46", "15" = "#7ab062",
  "16" = "#e27596",
  "6" = "#cb3e83", "10" = "#653f8c", "18" = "#9f5a2a", "22" = "#7c7cb0",
  "2" = "#95bcd9", "5" = "#2253a1", "8" = "#3579b5",
  "13" = "#814722", "20" = "#ef9f98",
  "17" = "#388045",
  "12" = "#cd6c9b", "24" = "#cd6c9b",
  "4" = "#479ebe", "21" = "#398287", "23" = "#80bcba",
  "19" = "#717071",
  "25" = "#c5a730"
)
cols2 <- c(
  "Keratinocyte" = "#cb8929",
  "Hair follicle" = "#7ab062",
  "Lymphoid" = "#018387",
  "Myeloid" = "#7f7db9",
  "Melanocyte" = "#727274",
  "Fibroblast" = "#368fbf",
  "Endothelial" = "#f49797",
  "Pericyte & Smooth muscle" = "#19964f",
  "Mast" = "#d866a6",
  "Schwann" = "#c5a730"
)
cols <- c(cols1, cols2)

# ---- 标签中位坐标 ----------------------------------------------------------
sub_type_med <- df %>%
  group_by(`RNA_snn_res.0.5`) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

# 大群名放在簇外空白，避免压住数字；不是原文 median+1。
main_type_med <- data.frame(
  newMainCellTypes = c(
    "Keratinocyte", "Hair follicle", "Schwann", "Pericyte & Smooth muscle",
    "Melanocyte", "Fibroblast", "Endothelial", "Myeloid", "Lymphoid", "Mast"
  ),
  umap_1 = c(1.2, -1.8, -7.8, -14.2, -5.6, -7.2, 0.8, 12.4, 12.8, 14.8),
  umap_2 = c(11.6, 11.4, 6.8, -1.2, -1.6, -11.4, -5.4, 2.2, -13.4, 7.6),
  stringsAsFactors = FALSE
)

# ---- 虚线非凸包 + 亚群点 + 标签 + 小箭头轴 --------------------------------
p_umap <- ggplot(df, aes(x = umap_1, y = umap_2)) +
  stat_unchull(
    aes(fill = newMainCellTypes, color = newMainCellTypes),
    alpha = 0.05, linewidth = 0.7, lty = 2, delta = 0.25, show.legend = FALSE
  ) +
  geom_point(aes(color = `RNA_snn_res.0.5`), size = 0.22, stroke = 0, show.legend = FALSE) +
  geom_text_repel(
    data = sub_type_med,
    aes(label = `RNA_snn_res.0.5`),
    color = "black",
    fontface = "bold",
    size = 3.6,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = main_type_med,
    aes(label = newMainCellTypes, color = newMainCellTypes),
    fontface = "bold",
    size = 5.2,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cols, guide = "none") +
  scale_fill_manual(values = cols, guide = "none") +
  coord_equal() +
  theme_dr(xlength = 0.22, ylength = 0.22) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(
      arrow = arrow(type = "closed", length = unit(0.08, "inches"))
    ),
    axis.title = element_text(hjust = 0.05, size = 11),
    plot.margin = margin(8, 28, 8, 8)
  )

# ---- 导出 ------------------------------------------------------------------
out_dir <- file.path(root, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(out_dir, "synthetic_render.png"),
  p_umap,
  width = 7.2,
  height = 7.2,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(out_dir, "synthetic_render.pdf"),
  p_umap,
  width = 7.2,
  height = 7.2,
  bg = "white"
)
file.copy(
  file.path(out_dir, "synthetic_render.png"),
  file.path(root, "preview.png"),
  overwrite = TRUE
)
cat("wrote", file.path(out_dir, "synthetic_render.png"), "\n")

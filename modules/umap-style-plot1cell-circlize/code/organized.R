# =============================================================================
# 环形 Circos UMAP（plot1cell 风格）
# organized 版本：中文分节导航 + 线性脚本，不封装函数
# =============================================================================
# original.R 依赖 plot1cell::plot_circlize。整理版用 ggplot 同心极坐标近似。
# 数据是合成肝脏细胞坐标，不是 02.UMAP.Rdata。
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) %in% c("code", "validation")) dirname(script_dir) else script_dir

plot_df <- read.csv(file.path(root, "data/example_umap.csv"), stringsAsFactors = FALSE)
cluster_colors <- c(
  "Neutrophil" = "#BCBD22",
  "B cell" = "#1F77B4",
  "NK cell" = "#17BECF",
  "T cells" = "#98DF8A",
  "Plasma" = "#AEC7E8",
  "Kupffer" = "#8C564B",
  "Macrophage" = "#7F7F7F",
  "Conventional Dendritic Cell" = "#2CA02C",
  "LSECs" = "#E377C2",
  "Hepatocytes" = "#9467BD",
  "Cholangiocyte" = "#FF7F0E",
  "Plasmacytoid dendritic cells" = "#FFBB78",
  "Hepatic stellate cell" = "#D62728"
)
group_colors <- c("Control" = "#1F77B4", "Treat" = "#FF7F0E")
plot_df$celltype <- factor(plot_df$celltype, levels = names(cluster_colors))

# ---- 把 UMAP 缩放到中心圆盘 ------------------------------------------------
cx <- mean(range(plot_df$UMAP_1))
cy <- mean(range(plot_df$UMAP_2))
span <- max(diff(range(plot_df$UMAP_1)), diff(range(plot_df$UMAP_2))) / 2
plot_df$nx <- (plot_df$UMAP_1 - cx) / span
plot_df$ny <- (plot_df$UMAP_2 - cy) / span
n_type <- nlevels(plot_df$celltype)
plot_df$theta_x <- 0.5 + (atan2(plot_df$ny, plot_df$nx) + pi) / (2 * pi) * n_type
plot_df$r <- 0.08 + 0.62 * pmin(sqrt(plot_df$nx^2 + plot_df$ny^2), 1)
label_df <- plot_df %>%
  group_by(celltype) %>%
  summarise(theta_x = median(theta_x), r = median(r), .groups = "drop")

ring <- plot_df %>%
  count(celltype, group, name = "n") %>%
  group_by(celltype) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()
cluster_ring <- ring %>%
  distinct(celltype) %>%
  mutate(
    xmin = as.numeric(celltype) - 0.5,
    xmax = as.numeric(celltype) + 0.5,
    ymin = 0.82,
    ymax = 0.94,
    fill_id = as.character(celltype)
  )
# 外圈按扇区切分 Control/Treat，而不是整圈单色。
group_ring <- ring %>%
  group_by(celltype) %>%
  arrange(group, .by_group = TRUE) %>%
  mutate(
    xmin = as.numeric(celltype) - 0.5 + cumsum(lag(pct, default = 0)),
    xmax = as.numeric(celltype) - 0.5 + cumsum(pct),
    ymin = 0.96,
    ymax = 1.10,
    fill_id = as.character(group)
  ) %>%
  ungroup()
rings <- bind_rows(
  cluster_ring[, c("xmin", "xmax", "ymin", "ymax", "fill_id")],
  group_ring[, c("xmin", "xmax", "ymin", "ymax", "fill_id")]
)
fill_values <- c(cluster_colors, group_colors)

p1 <- ggplot() +
  geom_rect(
    data = rings,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_id),
    color = "white", linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = theta_x, y = r, color = celltype),
    size = 0.22, alpha = 0.85, stroke = 0
  ) +
  geom_text(
    data = label_df,
    aes(x = theta_x, y = r, label = celltype),
    size = 2.3, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(values = fill_values, breaks = c("Control", "Treat"), name = "group") +
  scale_color_manual(values = cluster_colors, guide = "none") +
  coord_polar(start = -pi / n_type) +
  ylim(0, 1.12) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

out_dir <- file.path(root, "validation")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
png_path <- file.path(out_dir, "synthetic_render.png")
pdf_path <- file.path(out_dir, "synthetic_render.pdf")
preview_path <- file.path(root, "preview.png")
ggsave(filename = png_path, plot = p1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(filename = pdf_path, plot = p1, width = 8, height = 8, bg = "white")
file.copy(png_path, preview_path, overwrite = TRUE)
message("wrote ", png_path)

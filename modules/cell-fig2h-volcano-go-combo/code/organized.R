# =============================================================================
# Cell Fig. 2H 风格组合图：百分比差气泡散点 + 上下调 GO 条形图
# organized 版本：中文分节导航 + 线性脚本，不封装自写绘图函数
# =============================================================================
# 这不是论文或微信作者的原始绘图脚本。
# GitHub / Zenodo 没有 Figure2H；微信只公开了 PBMC 入口。
# 本脚本按微信文中的 Cell 原图截图重建绘图层，输入是合成表。
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

# ---- 第一步：定位条目根目录 ------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") {
  dirname(script_dir)
} else if (basename(script_dir) == "validation") {
  dirname(script_dir)
} else {
  script_dir
}

# ---- 第二步：读合成 DEG 与 GO 表 ------------------------------------------
de <- read.csv(file.path(root, "data/example_de.csv"), stringsAsFactors = FALSE)
go <- read.csv(file.path(root, "data/example_go.csv"), stringsAsFactors = FALSE)

col_up <- "#E9584B"
col_down <- "#3B86B6"
col_ns <- "#BDBDBD"
col_bar_up <- "#EF8079"
col_bar_down <- "#7EA7C8"
col_title <- "#0B4F86"

de$sig <- factor(de$sig, levels = c("NS", "Up", "Down"))
labeled <- subset(de, source == "labeled")

# ---- 第三步：左侧百分比差 × log2FC 气泡散点 --------------------------------
p_scatter <- ggplot(de, aes(x = pct_diff, y = avg_log2FC)) +
  geom_hline(yintercept = c(-1, 1), linetype = "longdash", linewidth = 0.62, colour = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "longdash", linewidth = 0.62, colour = "black") +
  geom_point(
    data = subset(de, source == "background"),
    aes(size = size_value),
    colour = col_ns,
    alpha = 0.72,
    stroke = 0
  ) +
  geom_point(
    data = labeled,
    aes(size = size_value, colour = sig),
    alpha = 0.96,
    stroke = 0
  ) +
  geom_segment(
    data = labeled,
    aes(
      x = pct_diff,
      y = avg_log2FC,
      xend = pct_diff + nudge_x * 0.58,
      yend = avg_log2FC + nudge_y * 0.58
    ),
    inherit.aes = FALSE,
    linewidth = 0.22,
    colour = "grey15"
  ) +
  geom_text(
    data = labeled,
    aes(x = pct_diff + nudge_x, y = avg_log2FC + nudge_y, label = gene),
    inherit.aes = FALSE,
    fontface = "italic",
    size = 2.35,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 11.35,
    label = "Co_ESCs & AMLCs\nVS\nESCs_only",
    colour = col_title,
    fontface = "bold",
    size = 4.35,
    lineheight = 0.90
  ) +
  scale_colour_manual(values = c(Up = col_up, Down = col_down), guide = "none") +
  scale_size_continuous(range = c(0.85, 4.8), guide = "none") +
  scale_x_continuous(
    limits = c(-0.84, 0.91),
    breaks = c(-0.5, 0.0, 0.5),
    labels = c("-0.5", "0.0", "0.5"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-6.1, 12.8),
    breaks = c(-5, 0, 5, 10),
    expand = c(0, 0)
  ) +
  labs(x = "Percentage Difference", y = "Log2-Fold Change") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.40),
    axis.ticks.length = unit(2.0, "pt"),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.62),
    plot.margin = margin(10, 4, 6, 16)
  )

# ---- 第四步：右侧颜色图例与点大小图例 --------------------------------------
legend_color <- data.frame(
  x = 0.08,
  y = c(0.58, 0.28),
  lab = c("Up regulation", "Down regulation"),
  col = c(col_up, col_down),
  stringsAsFactors = FALSE
)
legend_size <- data.frame(
  x = 1.28,
  y = seq(0.82, 0.10, length.out = 5),
  lab = c("0.0", "0.2", "0.4", "0.6", "0.8"),
  sz = c(1.2, 2.1, 2.9, 3.7, 4.5),
  stringsAsFactors = FALSE
)

p_legend <- ggplot() +
  geom_point(
    data = legend_color,
    aes(x = x, y = y),
    colour = legend_color$col,
    size = 3.7,
    stroke = 0
  ) +
  geom_text(
    data = legend_color,
    aes(x = x + 0.09, y = y, label = lab),
    hjust = 0,
    size = 3.15,
    colour = "black"
  ) +
  annotate(
    "text",
    x = 1.05,
    y = 0.99,
    label = "Percentage Difference",
    hjust = 0,
    vjust = 1,
    size = 3.15,
    colour = "black"
  ) +
  geom_point(
    data = legend_size,
    aes(x = x, y = y, size = sz),
    colour = "black",
    stroke = 0
  ) +
  geom_text(
    data = legend_size,
    aes(x = x + 0.14, y = y, label = lab),
    hjust = 0,
    size = 3.2,
    colour = "black"
  ) +
  scale_size_identity() +
  scale_x_continuous(limits = c(0, 2.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 6, 0, 0))

# ---- 第五步：右侧上下调 GO 条 ----------------------------------------------
# 原图条形从面板左缘起笔，文字叠在条上；坐标轴含 -10。
go$fill <- ifelse(go$direction == "Up", col_bar_up, col_bar_down)

p_go <- ggplot(go) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = y - 0.52, ymax = y + 0.52, fill = fill),
    colour = NA
  ) +
  geom_text(
    aes(x = xmin + 0.4, y = y, label = term),
    hjust = 0,
    size = 3.05,
    colour = "black"
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(-16, 32),
    breaks = c(-10, 0, 10, 20, 30),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(-1.0, 15.4), expand = c(0, 0)) +
  labs(x = "-log(P Value)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 0.40),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.40),
    plot.margin = margin(0, 8, 6, 0)
  )

# ---- 第六步：拼合并导出 ----------------------------------------------------
p_right <- p_legend / p_go + plot_layout(heights = c(0.30, 1))
p_fig <- (
  p_scatter + p_right + plot_layout(widths = c(1.22, 0.88))
) + plot_annotation(
  title = "H",
  theme = theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0,
      margin = margin(0, 0, -4, 0)
    ),
    plot.margin = margin(4, 8, 4, 4)
  )
)

out_dir <- file.path(root, "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
png_path <- file.path(out_dir, "synthetic_render.png")
pdf_path <- file.path(out_dir, "synthetic_render.pdf")
preview_path <- file.path(root, "preview.png")

ggsave(png_path, plot = p_fig, width = 7.7, height = 6.4, dpi = 300, bg = "white")
ggsave(pdf_path, plot = p_fig, width = 7.7, height = 6.4, bg = "white")
ggsave(preview_path, plot = p_fig, width = 7.7, height = 6.4, dpi = 220, bg = "white")

cat("wrote", png_path, "\n")
cat("wrote", pdf_path, "\n")
cat("wrote", preview_path, "\n")

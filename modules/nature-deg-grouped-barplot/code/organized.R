# =============================================================================
# Nature 风格分组 DEG 条形图（三时间点并排）
# organized：线性脚本 + 中文分节
# =============================================================================
# 与 original.R 的可见差异：
#   1. 读本草稿 data/DEGs.csv，输出写到 output/figures/
#   2. 三列间距从 20 改为 80，避免 “500” 与下一列 “0” 粘成 “5000”
#   3. 左边距加大，避免 Endothelial 被裁
#   4. PNG 用 ragg，PDF 用 cairo
# 数据：从微信文章数据表抄录的 23 个细胞类型 × 3 个 SNI 时间点 DEG 计数。
# 不是原文单细胞分析复现。
# =============================================================================

library(ggplot2)
library(readr)
library(ragg)

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
draft_dir <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_font_family <- "Arial"

# ---- 第一步：读 DEG 计数表 -------------------------------------------------
degs <- read_csv(file.path(draft_dir, "data", "DEGs.csv"), show_col_types = FALSE)
degs$cell_type <- factor(degs$cell_type, levels = rev(degs$cell_type))
n_type <- length(degs$cell_type)

Group_color <- c(
  "Intratelencephalic" = "#0072BD",
  "Pyramidal/Corticothalamic" = "#9E005D",
  "GABA-ergic" = "#00A69C",
  "Non-neuronal" = "#A3784D"
)
time_fill <- c(day = "#FF69B4", week = "#F72F8C", month = "#9F044D")

# 三列起点。原文 0/520/1040 会让 500 与 0 贴在一起。
col_w <- 500
gap <- 80
x0 <- c(0, col_w + gap, 2 * (col_w + gap))
names(x0) <- c("day", "week", "month")
x_right <- x0["month"] + col_w

n_nn <- sum(degs$Group == "Non-neuronal")
n_ga <- sum(degs$Group == "GABA-ergic")
n_py <- sum(degs$Group == "Pyramidal/Corticothalamic")

# ---- 第二步：三列条 + 分组底色 + 轴 ----------------------------------------
p <- ggplot(degs, aes(y = cell_type, colour = Group)) +
  labs(x = NULL, y = NULL) +
  geom_rect(
    aes(xmin = x0["day"], xmax = x0["day"] + SNI_3day,
        ymin = as.numeric(cell_type) - 0.3, ymax = as.numeric(cell_type) + 0.3),
    linewidth = 0, color = NA, fill = time_fill["day"]
  ) +
  geom_text(aes(x = x0["day"] - 12, label = SNI_3day, colour = Group),
            hjust = 1, size = 3.2, family = plot_font_family) +
  geom_rect(
    aes(xmin = x0["week"], xmax = x0["week"] + SNI_3week,
        ymin = as.numeric(cell_type) - 0.3, ymax = as.numeric(cell_type) + 0.3),
    linewidth = 0, color = NA, fill = time_fill["week"]
  ) +
  geom_text(aes(x = x0["week"] - 12, label = SNI_3week, colour = Group),
            hjust = 1, size = 3.2, family = plot_font_family) +
  geom_rect(
    aes(xmin = x0["month"], xmax = x0["month"] + SNI_3month,
        ymin = as.numeric(cell_type) - 0.3, ymax = as.numeric(cell_type) + 0.3),
    linewidth = 0, color = NA, fill = time_fill["month"]
  ) +
  geom_text(aes(x = x0["month"] - 12, label = SNI_3month, colour = Group),
            hjust = 1, size = 3.2, family = plot_font_family) +
  geom_segment(
    data = data.frame(x = unname(x0)),
    aes(x = x, xend = x, y = 0.5, yend = n_type + 0.5),
    inherit.aes = FALSE, linewidth = 0.55
  ) +
  geom_segment(
    data = data.frame(x = unname(x0)),
    aes(x = x, xend = x + col_w, y = n_type + 0.5, yend = n_type + 0.5),
    inherit.aes = FALSE, linewidth = 0.7
  ) +
  geom_segment(
    data = data.frame(x = c(unname(x0), unname(x0) + col_w)),
    aes(x = x, xend = x, y = n_type + 0.5, yend = n_type + 0.8),
    inherit.aes = FALSE, linewidth = 0.7
  ) +
  geom_text(
    data = data.frame(x = unname(x0)),
    aes(x = x, y = n_type + 1.05, label = "0"),
    size = 3.6, vjust = 0, inherit.aes = FALSE, family = plot_font_family
  ) +
  geom_text(
    data = data.frame(x = unname(x0) + col_w),
    aes(x = x, y = n_type + 1.05, label = "500"),
    size = 3.6, vjust = 0, hjust = 1, inherit.aes = FALSE, family = plot_font_family
  ) +
  annotate(
    "text", x = mean(c(0, x_right)), y = n_type + 2.15,
    label = "# Differentially Expressed Genes (DEGs)",
    hjust = 0.5, vjust = 0, size = 4.2, family = plot_font_family
  ) +
  annotate("segment", x = -110, xend = -110, y = 0.5, yend = n_type + 0.5, linewidth = 0.9) +
  geom_text(
    aes(x = -132, label = cell_type, colour = Group),
    hjust = 1, fontface = "bold", size = 3.3, family = plot_font_family
  ) +
  annotate("rect", xmin = -110, xmax = x_right,
           ymin = 0.5, ymax = n_nn + 0.45,
           alpha = 0.07, fill = unname(Group_color["Non-neuronal"])) +
  annotate("rect", xmin = -110, xmax = x_right,
           ymin = n_nn + 0.55, ymax = n_nn + n_ga + 0.45,
           alpha = 0.07, fill = unname(Group_color["GABA-ergic"])) +
  annotate("rect", xmin = -110, xmax = x_right,
           ymin = n_nn + n_ga + 0.55, ymax = n_nn + n_ga + n_py + 0.45,
           alpha = 0.07, fill = unname(Group_color["Pyramidal/Corticothalamic"])) +
  annotate("rect", xmin = -110, xmax = x_right,
           ymin = n_nn + n_ga + n_py + 0.55, ymax = n_type + 0.5,
           alpha = 0.07, fill = unname(Group_color["Intratelencephalic"])) +
  annotate("point", x = x0["day"], y = -1.15, shape = 15, size = 4.2, color = time_fill["day"]) +
  annotate("text", x = x0["day"] + 28, y = -1.15, label = "SNI 3-day",
           hjust = 0, vjust = 0.5, size = 3.6, family = plot_font_family) +
  annotate("point", x = x0["week"], y = -1.15, shape = 15, size = 4.2, color = time_fill["week"]) +
  annotate("text", x = x0["week"] + 28, y = -1.15, label = "SNI 3-week",
           hjust = 0, vjust = 0.5, size = 3.6, family = plot_font_family) +
  annotate("point", x = x0["month"], y = -1.15, shape = 15, size = 4.2, color = time_fill["month"]) +
  annotate("text", x = x0["month"] + 28, y = -1.15, label = "SNI 3-month",
           hjust = 0, vjust = 0.5, size = 3.6, family = plot_font_family) +
  scale_color_manual(values = Group_color) +
  theme_void(base_family = plot_font_family) +
  coord_cartesian(xlim = c(-430, x_right + 30), ylim = c(-2.4, n_type + 3.2), clip = "off") +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 12, r = 18, b = 16, l = 8)
  )

# ---- 第三步：导出 -----------------------------------------------------------
ggsave(file.path(out_dir, "nature_deg_grouped_barplot.png"), p,
       width = 8.6, height = 7.4, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "nature_deg_grouped_barplot.pdf"), p,
       width = 8.6, height = 7.4, device = cairo_pdf, bg = "white")
message("wrote ", file.path(out_dir, "nature_deg_grouped_barplot.png"))

# =============================================================================
# Nature Cell Biology Fig. 2d：HOX 山脊图 + 平均表达热图
# organized 版本：线性脚本 + 中文分节导航
# =============================================================================
# 与 code/original.R（微信原文）的可见差异：
#   1. 输入改为本草稿 data/ 下的官方 Source Data Fig. 2（MOESM7，sheet Panel d），
#      不再使用原文不存在的 references/...xlsx 路径。
#   2. 输出写到本草稿 output/figures/，PNG 用 ragg。
#   3. HOXB 热图色阶按期刊原图用 −3 到 3；原文四个簇一律 −4 到 4。
#   4. 其余山脊平滑、Hermite 尾、配色和 2×2 拼图保留原文逻辑。
# 数据来源：Tan et al., Nat Cell Biol 2026, Source Data Fig. 2
# =============================================================================

library(ggplot2)
library(ggridges)
library(patchwork)
library(cowplot)
library(scales)
library(ragg)

# ---- 路径：以本脚本所在目录为基准 ------------------------------------------
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      normalizePath(getwd())
    }
  }
)
if (basename(script_dir) == "code") {
  draft_dir <- dirname(script_dir)
} else {
  draft_dir <- script_dir
}

input_csv <- file.path(draft_dir, "data", "panel-d.csv")
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_png <- file.path(out_dir, "fig2d_hox_ridge_heatmap.png")
output_pdf <- file.path(out_dir, "fig2d_hox_ridge_heatmap.pdf")

plot_font_family <- "Arial"
# 比原文 7 in 略宽，给斜体基因名和两列中间留边。
figure_width_in <- 8.2
figure_height_in <- figure_width_in * 923 / 1438
figure_dpi <- 300
ridge_x_max <- 118
ridge_smooth_spar <- 0.45
ridge_tail_end <- 112

# ---- 第一步：读取官方 Source Data Panel d -----------------------------------
# 31 个 HOX 基因 × 100 个 A–P 空间层（AP_1–AP_100）。
# 数值是每层平均表达，不是单细胞 ridgeline 的原始观测。
panel_d <- as.data.frame(read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE))
gene_col <- names(panel_d)[1]
genes <- as.character(panel_d[[1]])
ap_names <- names(panel_d)[-1]
ap <- seq_along(ap_names)
expr <- as.matrix(panel_d[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- genes
colnames(expr) <- ap_names

groups <- list(
  HOXA = grep("^HOXA", genes),
  HOXB = grep("^HOXB", genes),
  HOXC = grep("^HOXC", genes),
  HOXD = grep("^HOXD", genes)
)

ridge_palettes <- list(
  HOXA = c("#c95665", "#b18fa6", "#c66e60", "#d0a16e", "#888b90",
           "#7e98bb", "#90aaa5", "#d09d68", "#8fa38d", "#eb5d60"),
  HOXB = c("#af91a9", "#b56d60", "#d2a570", "#888e92", "#7895b8",
           "#91a59a", "#d19c6c", "#8da382", "#eb5d60"),
  HOXC = c("#8b9091", "#7894b6", "#91aaa3", "#d7945d", "#8da17e", "#ed5b5e"),
  HOXD = c("#8b9091", "#7894b6", "#91aaa3", "#d7945d", "#8da17e", "#ed5b5e")
)

# 期刊原图 HOXB 色条是 −3/0/3，其余簇是 −4/0/4。
heat_limits <- list(
  HOXA = c(-4, 4),
  HOXB = c(-3, 3),
  HOXC = c(-4, 4),
  HOXD = c(-4, 4)
)

# ---- 第二步：行 z-score，只用于右侧热图 ------------------------------------
# 山脊高度另做 min-max，避免高表达基因压扁其他曲线。
row_z <- function(x) {
  z <- t(apply(x, 1, function(v) {
    s <- sd(v)
    if (is.na(s) || s == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  }))
  rownames(z) <- rownames(x)
  colnames(z) <- colnames(x)
  z
}

# ---- 第三步：单个 HOX 簇 = 山脊 + 热图 --------------------------------------
# 原文已有 make_panel()，整理版原样保留，不拆成更多绘图封装。
make_panel <- function(group_name, idx) {
  group_genes <- genes[idx]
  group_expr <- expr[idx, , drop = FALSE]
  z_lim <- heat_limits[[group_name]]

  ridge_height <- t(apply(group_expr, 1, function(v) {
    v <- v - min(v)
    if (max(v) == 0) return(rep(0, length(v)))
    v / max(v)
  }))
  rownames(ridge_height) <- group_genes
  colnames(ridge_height) <- ap_names

  ridge_df <- do.call(rbind, lapply(seq_along(group_genes), function(i) {
    # 100 层曲线轻度平滑；100 之后的尾巴只是显示延伸，不是新的 AP 测量。
    smooth_fit <- smooth.spline(
      ap, as.numeric(ridge_height[i, ]), spar = ridge_smooth_spar
    )
    smooth_x <- seq(min(ap), ridge_x_max, length.out = 360)
    smooth_y <- predict(smooth_fit, x = pmin(smooth_x, max(ap)))$y
    tail_idx <- smooth_x > max(ap)
    if (any(tail_idx)) {
      x0 <- max(ap)
      x1 <- ridge_tail_end
      y0 <- predict(smooth_fit, x = x0)$y
      y1 <- 0
      slope0 <- (predict(smooth_fit, x = x0)$y -
        predict(smooth_fit, x = x0 - 5)$y) / 5
      tail_x <- pmin(smooth_x[tail_idx], x1)
      tail_fraction <- (tail_x - x0) / (x1 - x0)
      tail_length <- x1 - x0
      smooth_y[tail_idx] <-
        (2 * tail_fraction^3 - 3 * tail_fraction^2 + 1) * y0 +
        (tail_fraction^3 - 2 * tail_fraction^2 + tail_fraction) *
          tail_length * slope0 +
        (-2 * tail_fraction^3 + 3 * tail_fraction^2) * y1
      smooth_y[smooth_x > x1] <- 0
    }
    data.frame(
      AP = smooth_x,
      gene = factor(group_genes[i], levels = rev(group_genes)),
      y = length(group_genes) - i + 1,
      height = pmax(0, smooth_y)
    )
  }))

  ridge <- ggplot(ridge_df, aes(x = AP, y = y, height = height, fill = gene)) +
    geom_ridgeline(
      scale = 1.15, min_height = 0, color = "#1f1f1f", linewidth = 0.18,
      alpha = 1
    ) +
    scale_fill_manual(
      values = setNames(
        ridge_palettes[[group_name]][seq_along(group_genes)],
        group_genes
      )
    ) +
    scale_x_continuous(
      breaks = c(1, 21, 41, 61, 81, 100),
      labels = c("0", "20", "40", "60", "80", "100"),
      limits = c(0, ridge_x_max),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(group_genes),
      labels = rev(group_genes),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    coord_cartesian(xlim = c(0, ridge_x_max), expand = FALSE, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 6, base_family = plot_font_family) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "italic", size = 6.2),
      axis.text.x = element_text(size = 5.2),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "#e5e5e5", linewidth = 0.15),
      plot.margin = margin(2, 6, 2, 10)
    )

  heat_z <- row_z(group_expr)
  heat_df <- do.call(rbind, lapply(seq_along(group_genes), function(i) {
    data.frame(
      gene = factor(group_genes[i], levels = rev(group_genes)),
      AP = ap,
      z = as.numeric(heat_z[i, ])
    )
  }))

  heat <- ggplot(heat_df, aes(x = AP, y = gene, fill = z)) +
    geom_tile(width = 1, height = 1) +
    scale_x_continuous(
      breaks = c(1, 100), labels = c("A", "P"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(position = "right", drop = FALSE) +
    scale_fill_gradient2(
      low = "#425794", mid = "#f3f3f3", high = "#e4515a", midpoint = 0,
      limits = z_lim, oob = squish,
      breaks = c(z_lim[1], 0, z_lim[2]),
      name = "Exp."
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 6, base_family = plot_font_family) +
    theme(
      axis.text.y = element_text(
        face = "italic", size = 6.2, hjust = 0, margin = margin(l = 2)
      ),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 5.2),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 4.8),
      legend.key.width = unit(0.35, "cm"),
      legend.key.height = unit(0.10, "cm"),
      legend.box.spacing = unit(0, "pt"),
      plot.margin = margin(2, 16, 2, 4)
    ) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(0.85, "cm"), barheight = unit(0.10, "cm"),
      ticks = FALSE
    ))

  ridge_title <- ggdraw() + draw_label(
    paste0(group_name, " clusters"), x = 0.5, y = 0.5,
    hjust = 0.5, vjust = 0.5, size = 9, fontfamily = plot_font_family
  )
  heat_title <- ggdraw() + draw_label(
    "Average expression level", x = 0.5, y = 0.5,
    hjust = 0.5, vjust = 0.5, size = 7, fontfamily = plot_font_family
  )
  heat_legend <- get_legend(heat)
  heat_no_legend <- heat + theme(legend.position = "none")
  panel_widths <- c(1.18, 0.14, 1.00)
  blank <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  title_row <- plot_grid(
    ridge_title, blank, heat_title,
    ncol = 3, rel_widths = panel_widths
  )
  body_row <- plot_grid(
    ridge, blank, heat_no_legend,
    ncol = 3, rel_widths = panel_widths
  )
  ridge_footer <- ggdraw() + draw_label(
    "A  ->  P", x = 0.5, y = 0.5,
    hjust = 0.5, vjust = 0.5, size = 6.0, fontfamily = plot_font_family
  )
  heat_ap_footer <- ggdraw() + draw_label(
    "A  ->  P", x = 0.38, y = 0.5,
    hjust = 0.5, vjust = 0.5, size = 6.0, fontfamily = plot_font_family
  )
  heat_footer <- plot_grid(
    heat_ap_footer, heat_legend,
    ncol = 2, rel_widths = c(0.78, 1.22)
  )
  footer_row <- plot_grid(
    ridge_footer, blank, heat_footer,
    ncol = 3, rel_widths = panel_widths
  )
  plot_grid(
    title_row, body_row, footer_row,
    ncol = 1, rel_heights = c(0.15, 1, 0.16), align = "v"
  )
}

# ---- 第四步：2×2 拼图并导出 ------------------------------------------------
# 以上完成四个簇的单独面板；接下来按期刊 Fig. 2d 排成 HOXA/B 上、HOXC/D 下。
fig <- wrap_plots(
  make_panel("HOXA", groups$HOXA),
  make_panel("HOXB", groups$HOXB),
  make_panel("HOXC", groups$HOXC),
  make_panel("HOXD", groups$HOXD),
  ncol = 2
) +
  plot_annotation(
    title = "d",
    theme = theme(
      plot.title = element_text(
        hjust = 0, face = "bold", size = 11, family = plot_font_family
      ),
      plot.margin = margin(6, 14, 8, 10)
    )
  )

ggsave(
  filename = output_png, plot = fig,
  width = figure_width_in, height = figure_height_in,
  dpi = figure_dpi, bg = "white", device = ragg::agg_png
)
ggsave(
  filename = output_pdf, plot = fig,
  width = figure_width_in, height = figure_height_in,
  device = cairo_pdf, bg = "white"
)

message("wrote ", output_png)
message("wrote ", output_pdf)

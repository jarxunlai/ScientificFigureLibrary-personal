# =============================================================================
# Science 风格 Mantel 网络相关热图
# organized：作者微信代码的可运行版（linkET::qcorrplot + geom_couple）
# =============================================================================
# 与 code/original.R（微信原文）的可见差异：
#   1. 从 vegan 写出 data/varespec.csv、data/varechem.csv（作者本地
#      ./dataset/ 未随文提供；列数正好是 44 与 14，对得上 spec_select）。
#   2. 输出写到本草稿 output/figures/，PNG 用 ragg。
#   3. 删掉原文粘进来的 SHAP 段和 p2 气泡图（那是另一篇 Python 文）。
#   4. 绘图函数保持作者 p1：type = "lower", diag = FALSE, geom_square,
#      geom_mark, geom_couple。这是 linkET 原样，不再用 ggplot 手动画三角。
#   5. 原文 library(tidyverse)；本环境未装元包，改为 dplyr（%>% / mutate）。
# =============================================================================

library(linkET)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggnewscale)
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
data_dir <- file.path(draft_dir, "data")
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 第一步：读组成表与环境表 ----------------------------------------------
# vegan::varespec 24×44 地衣物种；vegan::varechem 24×14 土壤化学。
# 作者用这套表演示 Mantel 连线，不是 Tara Oceans 原表。
varespec <- read.csv(file.path(data_dir, "varespec.csv"))
varechem <- read.csv(file.path(data_dir, "varechem.csv"))

# ---- 第二步：Mantel 检验并切成图例箱子 --------------------------------------
mantel <- mantel_test(
  varespec, varechem,
  spec_select = list(
    `Taxonomic\ncomposition\n(16S OTUs)` = 1:18,
    `Gene\nfunctional\ncomposition` = 19:32,
    `Taxonomic\ncomposition\n(mOTUs)` = 33:44
  )
) %>%
  mutate(
    rd = cut(r, breaks = c(-Inf, 0.2, 0.3, Inf),
             labels = c("< 0.2", "0.2 - 0.3", ">= 0.3")),
    pd = cut(p, breaks = c(-Inf, 0.005, 0.01, 0.05, Inf),
             labels = c("< 0.005", "0.005 - 0.01", "0.01 - 0.05", ">= 0.05"))
  )

# ---- 第三步：作者 p1，qcorrplot + geom_couple -------------------------------
# type = "lower" 按原文。linkET 的 lower 在成品里常表现为对角左下填充、
# 组成节点在对角另一侧；这是包函数布局，不再手动画上/下三角。
p1 <- qcorrplot(
  correlate(varechem),
  grid_col = "grey50",
  grid_size = 0.5,
  type = "lower",
  diag = FALSE
) +
  geom_square() +
  geom_mark(
    size = 4,
    only_mark = TRUE,
    sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05,
    colour = "white"
  ) +
  geom_couple(
    data = mantel,
    aes(color = pd, size = rd),
    label.size = 3.88,
    label.family = "",
    label.fontface = 1,
    nudge_x = 0.2,
    curvature = nice_curvature(by = "from"),
    point_fill = "white",
    point_color = "gray50"
  ) +
  scale_fill_gradient2(
    limits = c(-0.8, 0.8),
    mid = "white",
    low = "#2a7bbb",
    high = "#e67e22",
    breaks = seq(-0.8, 0.8, 0.4)
  ) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_color_manual(values = c("#2d5d7a", "#5a90b3", "#f39c12", "#cccccc")) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      order = 2,
      keyheight = unit(0.5, "cm")
    ),
    colour = guide_legend(
      title = "Mantel's p",
      order = 1,
      keyheight = unit(0.5, "cm")
    ),
    fill = guide_colorbar(
      title = "Pearson's r",
      keyheight = unit(2.2, "cm"),
      keywidth = unit(0.5, "cm"),
      order = 3
    )
  ) +
  theme(legend.box.spacing = unit(0, "pt"))

# ---- 第四步：导出 -----------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "mantel_corrplot_square.png"),
  plot = p1, width = 8.8, height = 6, dpi = 300,
  bg = "white", device = ragg::agg_png
)
ggsave(
  filename = file.path(out_dir, "mantel_corrplot_square.pdf"),
  plot = p1, width = 8.8, height = 6,
  device = cairo_pdf, bg = "white"
)
message("wrote ", file.path(out_dir, "mantel_corrplot_square.png"))

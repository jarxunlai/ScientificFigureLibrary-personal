# =============================================================================
# GO 富集彗星图（叠加版）
# organized 版本：线性脚本 + 中文分节导航
# =============================================================================
# 与 code/original.R（微信原文）的可见差异：
#   1. 不重跑 GSE128033 / Seurat / clusterProfiler；从视觉重建长表读入；
#   2. 本条目只导出叠加图，不分面；
#   3. geom_text 用 aes(label = Count)；
#   4. linewidth 替代已弃用的 size。
# 数据来源：微信教程分面图读数重建，不是论文 Source Data，也不是 enrichGO 重算。
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(ggforce)
library(ragg)

# ---- 路径：以本脚本所在目录为基准 ------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir

input_csv <- file.path(root, "data", "enrichment_comet.csv")
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 第一步：读取视觉重建的富集长表 ----------------------------------------
dt <- read_csv(input_csv, show_col_types = FALSE)

cluster_levels <- c(
  "Macrophages",
  "ILC1/NK cells",
  "Endothelial cells",
  "Epithelial cells",
  "T cells",
  "Fibroblasts",
  "Mast cells",
  "Smooth muscle cells",
  "Cycling macrophages",
  "B cells",
  "Lymphatic endothelial cells"
)
dt <- dt %>%
  mutate(Cluster = factor(Cluster, levels = cluster_levels)) %>%
  arrange(pvalue, Cluster)
dt$Description <- factor(dt$Description, levels = unique(dt$Description))

# ---- 第二步：细胞类型配色（微信原文） --------------------------------------
colors <- c(
  "Macrophages" = "#B0C4DE",
  "Mast cells" = "#FF69B4",
  "Endothelial cells" = "#FFB6C1",
  "ILC1/NK cells" = "#FFFFE0",
  "B cells" = "#ADD8E6",
  "T cells" = "#90EE90",
  "Fibroblasts" = "#FFA07A",
  "Smooth muscle cells" = "#D3D3D3",
  "Epithelial cells" = "#D8BFD8",
  "Cycling macrophages" = "#E6E6FA",
  "Lymphatic endothelial cells" = "#FFD700"
)

# ---- 第三步：叠加彗星图 ----------------------------------------------------
p1 <- ggplot(dt) +
  geom_link(
    aes(
      x = 0,
      y = Description,
      xend = -log10(pvalue),
      yend = Description,
      alpha = after_stat(index),
      color = Cluster,
      linewidth = after_stat(index)
    ),
    n = 500,
    show.legend = TRUE
  ) +
  geom_point(
    aes(x = -log10(pvalue), y = Description),
    color = "black",
    fill = "white",
    size = 6,
    shape = 21
  ) +
  geom_text(
    aes(x = -log10(pvalue), y = Description, label = Count),
    size = 3,
    nudge_x = 0.05
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 13)
  ) +
  xlab("-Log10 Pvalue") +
  ylab("") +
  scale_color_manual(values = colors, breaks = cluster_levels)

# ---- 第四步：导出 ----------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "enrich_comet_combined.png"),
  plot = p1,
  width = 15,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  filename = file.path(out_dir, "enrich_comet_combined.pdf"),
  plot = p1,
  width = 15,
  height = 8,
  bg = "white"
)
ggsave(
  filename = file.path(root, "preview.png"),
  plot = p1,
  width = 15,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("wrote combined comet plot")

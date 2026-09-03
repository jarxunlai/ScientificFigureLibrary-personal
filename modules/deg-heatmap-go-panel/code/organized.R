# =============================================================================
# 差异基因热图 + GO 富集组合 Panel（Cell Discovery 2026 T1DM 复现）
# organized 版本：中文导航注释 + 规范缩进
# =============================================================================
# 与 code/original.R 的差异（仅此两类，代码逻辑零改动）：
#   1. 注释：来源/执行状态说明改为中文，并增加中文分节导航注释；
#   2. 排版：字符串等行的缩进规范化，便于逐行审阅。
# 任何影响输出的改动（排序、颜色、标签、参数默认值）均未做。
#
# 来源：微信公众号文章《高分Panel复现系列｜新玩法！把差异基因热图和 GO 富集
#   放进同一个 Panel》（公众号：多线程核糖体，2026-08），代码由用户从原文复制。
#   https://mp.weixin.qq.com/s/ASViS9MkzqbuQAlhaGudKQ
# 参考图：Single-cell atlas of reproductive endocrine organs reveals
#   transcriptomic responses to type 1 diabetes mellitus in nonhuman primates.
#   Cell Discovery, 2026（panel C）
#
# 执行状态：本项目未运行。文章未随附示例数据，需自备 3 个输入 CSV：
#   input_deg_heatmap.csv      列: gene, module, cell_type, tissue, regulation
#   input_cell_annotation.csv  列: tissue, cell_type, cell_index
#   input_go_enrichment.csv    列: module, term, neg_log10_p
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(grid)
library(scales)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE)) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir

# ---- 第一步：读取输入数据 --------------------------------------------------
# 三个 CSV 的列结构见文件头注释；更换数据时只需保证列名一致。
deg <- read_csv(
  file.path(root, "data", "input_deg_heatmap.csv"),
  show_col_types = FALSE
)

cell_anno <- read_csv(
  file.path(root, "data", "input_cell_annotation.csv"),
  show_col_types = FALSE
)

go <- read_csv(
  file.path(root, "data", "input_go_enrichment.csv"),
  show_col_types = FALSE
)

# ---- 第二步：固定基因与细胞类型顺序 ----------------------------------------
# 顺序必须提前固定，否则换数据重跑时布局会漂移；基因按模块（up/down/Other）分组排列。
cell_order <- cell_anno$cell_type

gene_order <- deg |>
  distinct(gene, module) |>
  mutate(
    module = factor(
      module,
      levels = c("Common up", "Common down", "Other")
    )
  ) |>
  arrange(module, gene) |>
  pull(gene)

plot_dat <- deg |>
  mutate(
    cell_type = factor(cell_type, levels = cell_order),
    gene = factor(gene, levels = rev(gene_order)),
    regulation = factor(
      regulation,
      levels = c(
        "Upregulated in T1DM",
        "Unchanged",
        "Downregulated in T1DM"
      )
    )
  )

# module_df：各模块在热图行方向上的起止位置，供左侧模块条使用。
module_df <- plot_dat |>
  distinct(gene, module) |>
  mutate(y = as.numeric(gene)) |>
  group_by(module) |>
  summarise(
    ymin = min(y) - 0.5,
    ymax = max(y) + 0.5,
    .groups = "drop"
  )

# ---- 第三步：配色 ------------------------------------------------------------
# reg_cols：差异状态配色（粉=上调，白=不变，蓝=下调）；tissue_cols：组织注释条独立配色，
# 两者故意分开，避免组织信息与差异状态混在同一套颜色里。
reg_cols <- c(
  "Upregulated in T1DM" = "#f23aa5",
  "Unchanged" = "#ffffff",
  "Downregulated in T1DM" = "#85c8e4"
)

tissue_cols <- c(
  "Hypothalamus" = "#f5cf3d",
  "Pituitary" = "#24a8df",
  "Ovary" = "#f07557",
  "Uterus" = "#70c644"
)

theme_clean <- theme_void(base_size = 8) +
  theme(plot.margin = margin(0, 0, 0, 0))

legend_dat <- tibble(
  x = c(0.02, 0.42, 0.68),
  label = names(reg_cols),
  fill = unname(reg_cols)
)

# ---- 第四步：自绘图例行 ------------------------------------------------------
# 不用 ggplot 默认 guide，而是用 geom_tile + geom_text 手工拼一行紧凑图例。
p_legend <- ggplot(
  legend_dat,
  aes(x = x, y = 0.52)
) +
  geom_tile(
    aes(fill = label),
    width = 0.034,
    height = 0.20,
    color = "grey40",
    linewidth = 0.13
  ) +
  geom_text(
    aes(
      x = x + 0.026,
      label = label
    ),
    hjust = 0,
    size = 1.85,
    color = "black"
  ) +
  scale_fill_manual(
    values = reg_cols,
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_clean

# p_tissue：顶部组织注释条，组织名标签居中放在各色块上。
p_tissue <- ggplot(
  cell_anno,
  aes(
    x = factor(cell_type, levels = cell_order),
    y = 1,
    fill = tissue
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.22
  ) +
  geom_text(
    data = cell_anno |>
      group_by(tissue) |>
      summarise(
        x = mean(
          as.numeric(
            factor(cell_type, levels = cell_order)
          )
        ),
        .groups = "drop"
      ),
    aes(
      x = x,
      y = 1,
      label = tissue
    ),
    inherit.aes = FALSE,
    size = 1.75,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(
    values = tissue_cols,
    guide = "none"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void(base_size = 8) +
  theme(plot.margin = margin(0, 0, 0, 0))

# ---- 第五步：主热图 ----------------------------------------------------------
# 两条白色 geom_hline 在 Common down / Common up 模块边界处留出分隔缝。
p_heat <- ggplot(
  plot_dat,
  aes(cell_type, gene, fill = regulation)
) +
  geom_tile(
    color = "#f2f2f2",
    linewidth = 0.02
  ) +
  geom_hline(
    yintercept = max(
      plot_dat |>
        filter(module == "Common down") |>
        mutate(y = as.numeric(gene)) |>
        pull(y)
    ) + 0.5,
    color = "white",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = min(
      plot_dat |>
        filter(module == "Common up") |>
        mutate(y = as.numeric(gene)) |>
        pull(y)
    ) - 0.5,
    color = "white",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = reg_cols,
    guide = "none"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 4.8,
      color = "black"
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- 第六步：左侧模块标记条 ---------------------------------------------------
# 注意：annotate 里的基因数（438/416）是作者示例数据的硬编码值，换数据后必须手改。
p_side <- ggplot(module_df) +
  geom_rect(
    aes(
      xmin = 0.28,
      xmax = 0.42,
      ymin = ymin,
      ymax = ymax,
      fill = module
    ),
    color = NA
  ) +
  geom_segment(
    aes(
      x = 0.50,
      xend = 0.50,
      y = ymin,
      yend = ymax,
      color = module
    ),
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = 0.15,
    y = mean(
      c(
        min(module_df$ymin[module_df$module == "Common up"]),
        max(module_df$ymax[module_df$module == "Common up"])
      )
    ),
    label = "438",
    color = "#f23aa5",
    size = 2.15,
    angle = 90
  ) +
  annotate(
    "text",
    x = 0.15,
    y = mean(
      c(
        min(module_df$ymin[module_df$module == "Common down"]),
        max(module_df$ymax[module_df$module == "Common down"])
      )
    ),
    label = "416",
    color = "#4fa7cc",
    size = 2.15,
    angle = 90
  ) +
  annotate(
    "text",
    x = 0.68,
    y = mean(
      c(
        min(module_df$ymin[module_df$module == "Common up"]),
        max(module_df$ymax[module_df$module == "Common up"])
      )
    ),
    label = "Common",
    color = "#f23aa5",
    size = 2.35,
    angle = 90
  ) +
  annotate(
    "text",
    x = 0.68,
    y = mean(
      c(
        min(module_df$ymin[module_df$module == "Common down"]),
        max(module_df$ymax[module_df$module == "Common down"])
      )
    ),
    label = "Common",
    color = "#4fa7cc",
    size = 2.35,
    angle = 90
  ) +
  scale_fill_manual(
    values = c(
      "Common up" = "#f23aa5",
      "Common down" = "#85c8e4",
      "Other" = "white"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "Common up" = "#f23aa5",
      "Common down" = "#4fa7cc",
      "Other" = "white"
    ),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = range(module_df$ymin, module_df$ymax),
    expand = FALSE
  ) +
  theme_void(base_size = 8)

# ---- 第七步：GO 富集条形图（上下两组共用一个绘图函数） ------------------------
# 关键手法：GO 条目文字不作为 y 轴标签，而是用 geom_text 直接叠加在条形内部。
go_plot <- function(df, fill_col) {
  df <- df |>
    mutate(
      term = as.character(term),
      y = rev(seq_along(term))
    )

  ggplot(df, aes(y = y)) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = neg_log10_p,
        ymin = y - 0.34,
        ymax = y + 0.34
      ),
      fill = fill_col,
      color = NA
    ) +
    geom_text(
      aes(
        x = 0.28,
        y = y,
        label = term
      ),
      hjust = 0,
      vjust = 0.5,
      size = 1.62,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(0, 25),
      breaks = seq(0, 25, 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0.45, max(df$y) + 0.55),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    labs(
      x = "-log10(Pvalue)",
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(
        size = 6.3,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 7.2,
        color = "black",
        margin = margin(t = 1)
      ),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.3),
      axis.ticks.x = element_line(linewidth = 0.3),
      plot.margin = margin(0, 0, 0, 0)
    )
}

p_go_up <- go_plot(
  go |> filter(module == "Common up"),
  "#f23aa5"
)

p_go_down <- go_plot(
  go |> filter(module == "Common down"),
  "#85c8e4"
)

p_go <- p_go_up / p_go_down

# ---- 第八步：热图与 GO 图之间的半透明连接带 ----------------------------------
# connector 生成楔形 polygonGrob，ytop/ybottom 是相对连接区域的比例位置（0–1）。
connector <- function(fill_col, ytop, ybottom) {
  polygonGrob(
    x = unit(c(0, 1, 1, 0), "npc"),
    y = unit(
      c(ytop, 0.78, 0.52, ybottom),
      "npc"
    ),
    gp = gpar(
      fill = alpha(fill_col, 0.24),
      col = NA
    )
  )
}

p_link <- ggplot() +
  annotation_custom(
    connector("#f23aa5", 0.93, 0.51),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  annotation_custom(
    connector("#85c8e4", 0.49, 0.07),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  theme_void()

# ---- 第九步：patchwork 组合并导出 --------------------------------------------
# 左侧块 = 图例行 / 组织条行 / 热图行；整体 = 左块 + 连接带 + GO 条形图，三栏宽比 1.08:0.15:0.78。
legend_row <- plot_spacer() + p_legend + plot_spacer() +
  plot_layout(widths = c(0.08, 1, 0.02))

tissue_row <- plot_spacer() + p_tissue + plot_spacer() +
  plot_layout(widths = c(0.08, 1, 0.02))

heatmap_row <- p_side + p_heat +
  plot_layout(widths = c(0.10, 1))

left_block <- legend_row / tissue_row / heatmap_row +
  plot_layout(
    heights = c(0.10, 0.07, 1)
  )

dir.create(file.path(root, "output"), recursive = TRUE, showWarnings = FALSE)

p <- wrap_plots(
  left_block,
  p_link,
  p_go,
  ncol = 3,
  widths = c(1.08, 0.15, 0.78)
)

ggsave(
  file.path(root, "preview.png"),
  p,
  width = 7.2,
  height = 4.35,
  dpi = 450,
  bg = "white"
)

ggsave(
  file.path(root, "output", "deg_heatmap_go_panel.pdf"),
  p,
  width = 7.2,
  height = 4.35,
  bg = "white"
)

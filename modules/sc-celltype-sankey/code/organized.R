# sc-celltype-sankey
#
# Clean, portable plotting example for Open Figure Modules. This linear R
# script reads only the synthetic CSV distributed with the module. External
# source images and unredistributed original code are not part of this module.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggforce)
})

WIDTH_IN <- 8
HEIGHT_IN <- 5

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
} else {
  script_dir
}

# ---- 第二步：读取合成细胞计数表 --------------------------------------------
counts <- read.csv(file.path(root, "data/example-counts-csv.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("Sample_ID", "Group", "celltype", "n") %in% names(counts)))
counts$Group <- factor(counts$Group, levels = c("Control", "Treatment"))

# ---- 第三步：细胞类型 × 分组流量 ------------------------------------------
cell_levels <- counts %>%
  group_by(celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(desc(n)) %>%
  pull(celltype)
links <- counts %>%
  group_by(celltype, Group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    celltype = factor(celltype, levels = cell_levels),
    Group = factor(Group, levels = c("Treatment", "Control"))
  )
links <- as.data.frame(links)
sets <- ggforce::gather_set_data(links, 1:2)

# 左侧节点 / 连线：细胞类型色。右侧节点：分组色。
pal_cell <- c(
  "Fibroblasts"         = "#85d2e1",
  "Endothelial cells"   = "#879da2",
  "smooth muscle cells" = "#02401b",
  "myeloid cells"       = "#e38501",
  "T-cells"             = "#e0d200",
  "Pericyte"            = "#eb5b5e",
  "NK-cells"            = "#ead2bd",
  "B-cells"             = "#749f87",
  "cardiomyocyte cells" = "#0073c2"
)
pal_group <- c(
  "Treatment"      = "#B2182B",
  "Control" = "#7092d2"
)
pal_all <- c(pal_cell, pal_group)

# ---- 第四步：先取出节点坐标，再画桑基 --------------------------------------
# ggforce 默认轴填充是深灰，这里改成 fill = y，并单独标注左右标签。
axis_width <- 0.13
x_expand <- expansion(add = c(0.95, 0.70))
p_axes <- ggplot(sets, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets_axes(axis.width = axis_width) +
  scale_x_discrete(expand = x_expand)
axis_df <- ggplot_build(p_axes)$data[[1]]
nodes <- axis_df %>%
  distinct(label, xmin, xmax, ymin, ymax) %>%
  mutate(
    ymid = (ymin + ymax) / 2,
    xmid = (xmin + xmax) / 2,
    is_left = xmid < 1.5,
    x_lab = ifelse(is_left, xmin - 0.04, xmax + 0.05),
    hjust = ifelse(is_left, 1, 0)
  )
group_key <- data.frame(
  Group = factor(c("Treatment", "Control"), levels = c("Treatment", "Control"))
)

# ---- 第五步：连线 + 着色节点 + 分组图例 ------------------------------------
p1 <- ggplot(sets, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(
    aes(fill = celltype),
    alpha = 0.55,
    colour = "grey25",
    linewidth = 0.25,
    axis.width = axis_width,
    strength = 0.5
  ) +
  geom_parallel_sets_axes(
    aes(fill = y),
    colour = "grey15",
    linewidth = 0.35,
    axis.width = axis_width
  ) +
  geom_text(
    data = nodes,
    aes(x = x_lab, y = ymid, label = label, hjust = hjust),
    inherit.aes = FALSE,
    size = 3,
    colour = "grey15"
  ) +
  geom_point(
    data = group_key,
    aes(x = 1, y = 0, colour = Group),
    inherit.aes = FALSE,
    shape = 15,
    size = 0,
    alpha = 0
  ) +
  scale_fill_manual(
    values = pal_all,
    name = "celltype",
    breaks = cell_levels
  ) +
  scale_colour_manual(
    values = pal_group,
    name = "group"
  ) +
  scale_x_discrete(
    labels = c("celltype", "group"),
    expand = x_expand
  ) +
  guides(
    fill = guide_legend(order = 1),
    colour = guide_legend(
      order = 2,
      override.aes = list(size = 5, alpha = 1, shape = 15)
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ---- 最后一步：导出生成预览 ------------------------------------------------
output_args <- sub("^--output-dir=", "", args[grepl("^--output-dir=", args)])
output_env <- Sys.getenv("SFL_OUTPUT_DIR", unset = "")
out_dir <- if (length(output_args) && nzchar(output_args[[1]])) {
  output_args[[1]]
} else if (nzchar(output_env)) {
  output_env
} else {
  file.path(tempdir(), "sc-celltype-sankey-output")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
preview_output <- Sys.getenv("SFL_PREVIEW_OUTPUT", unset = "")
render_path <- file.path(out_dir, "render.png")
ggsave(filename = render_path, plot = p1, width = WIDTH_IN, height = HEIGHT_IN, dpi = 300, bg = "white")
if (nzchar(preview_output)) {
  dir.create(dirname(preview_output), recursive = TRUE, showWarnings = FALSE)
  file.copy(render_path, preview_output, overwrite = TRUE)
}
message("wrote ", render_path)

# sc-celltype-stacked-area
#
# Clean, portable plotting example for Open Figure Modules. This linear R
# script reads only the synthetic CSV distributed with the module. External
# source images and unredistributed original code are not part of this module.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(scales)
})

WIDTH_IN <- 11
HEIGHT_IN <- 6.5

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
pal_npg <- c(
  "Fibroblasts"         = "#E64B35",
  "Endothelial cells"   = "#4DBBD5",
  "smooth muscle cells" = "#00A087",
  "myeloid cells"       = "#3C5488",
  "T-cells"             = "#F39B7F",
  "Pericyte"            = "#8491B4",
  "NK-cells"            = "#91D1C2",
  "B-cells"             = "#DC0000",
  "cardiomyocyte cells" = "#7E6148"
)

# ---- 第三步：样本顺序与比例 ----------------------------------------------
sample_levels <- counts %>%
  distinct(Sample_ID, Group) %>%
  arrange(Group, Sample_ID) %>%
  pull(Sample_ID)
prop_df <- counts %>%
  mutate(
    Sample_ID = factor(Sample_ID, levels = sample_levels),
    celltype = factor(celltype, levels = names(pal_npg))
  ) %>%
  group_by(Sample_ID) %>%
  mutate(pct = n / sum(n), idx = as.numeric(Sample_ID)) %>%
  ungroup()

# ---- 第四步：堆叠面积 + Control/HF 分隔线 --------------------------------
p1 <- ggplot(prop_df, aes(x = idx, y = pct, fill = celltype, group = celltype)) +
  # ggplot2 3.5 默认第一水平在顶部；reverse=TRUE 让 pal_npg 的 Fibroblasts 回到底部。
  geom_area(position = position_fill(reverse = TRUE), color = "white", linewidth = 0.3, alpha = 0.92) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  scale_x_continuous(
    breaks = seq_along(sample_levels),
    labels = sample_levels,
    expand = expansion(add = c(0.5, 0.8))
  ) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_manual(values = pal_npg, name = "Cell type",
                    guide = guide_legend(reverse = TRUE)) +
  labs(
    title = "Cell Type Composition across Samples",
    subtitle = "Stacked area chart (Control -> Treatment)",
    x = "Sample ID",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "grey30", face = "bold"),
    axis.text.y = element_text(color = "grey30"),
    axis.title = element_text(face = "bold", color = "grey20"),
    plot.title = element_text(face = "bold", size = 16, color = "grey15"),
    plot.subtitle = element_text(color = "grey50", size = 12),
    legend.title = element_text(face = "bold", color = "grey20"),
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
  file.path(tempdir(), "sc-celltype-stacked-area-output")
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

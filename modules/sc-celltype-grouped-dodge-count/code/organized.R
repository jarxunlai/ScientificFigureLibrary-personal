# sc-celltype-grouped-dodge-count
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

WIDTH_IN <- 10
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

# ---- 第三步：组 × 细胞类型绝对数 ------------------------------------------
count_df <- counts %>%
  group_by(Group, celltype) %>%
  summarise(n = sum(n), .groups = "drop")
ord <- counts %>%
  group_by(celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(n) %>%
  pull(celltype)
count_df$celltype <- factor(count_df$celltype, levels = ord)
pal_group <- c("Control" = "#2166AC", "Treatment" = "#B2182B")
dodge_width <- 0.7

# ---- 第四步：分组 dodge 柱 ------------------------------------------------
p1 <- ggplot(count_df, aes(x = celltype, y = n, fill = Group)) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = dodge_width,
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = comma(n)),
    position = position_dodge(width = dodge_width),
    vjust = -0.35,
    size = 3.3,
    color = "grey25",
    fontface = "bold"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = pal_group, name = "Group") +
  labs(
    title = "Absolute Cell Counts by Cell Type",
    subtitle = "Dodged by group (Control vs Treatment)",
    x = "Cell type",
    y = "Cell count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "grey30"),
    axis.title = element_text(face = "bold", color = "grey20"),
    plot.title = element_text(face = "bold", size = 16, color = "grey15"),
    plot.subtitle = element_text(color = "grey50", size = 12),
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
  file.path(tempdir(), "sc-celltype-grouped-dodge-count-output")
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

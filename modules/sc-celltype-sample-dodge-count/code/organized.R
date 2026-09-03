# sc-celltype-sample-dodge-count
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

WIDTH_IN <- 12
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

# ---- 第三步：样本 × 细胞类型绝对数 ----------------------------------------
count_df <- counts %>%
  group_by(Group, Sample_ID, celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(total_n = sum(n)) %>%
  ungroup()
ord <- counts %>%
  group_by(celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(n) %>%
  pull(celltype)
count_df$celltype <- factor(count_df$celltype, levels = ord)
sample_order <- count_df %>%
  distinct(Sample_ID, Group, total_n) %>%
  arrange(Group, desc(total_n)) %>%
  pull(Sample_ID)
count_df$Sample_ID <- factor(count_df$Sample_ID, levels = sample_order)
pal_samples <- c(
  "C1" = "#2166AC",
  "C2" = "#67A9CF",
  "T1" = "#67001F",
  "T2" = "#B2182B",
  "T3" = "#D6604D",
  "T4" = "#F4A582",
  "T5" = "#FDDBC7"
)
dodge_width <- 0.8

# ---- 第四步：分面 dodge 柱 ------------------------------------------------
p1 <- ggplot(count_df, aes(x = celltype, y = n, fill = Sample_ID)) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = dodge_width, color = "white", linewidth = 0.25
  ) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = pal_samples, name = "Sample") +
  labs(
    title = "Absolute Cell Counts by Cell Type",
    subtitle = "Dodged by sample across Control and Treatment groups",
    x = "Cell type",
    y = "Cell count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.spacing.x = unit(0.8, "cm"),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6),
    strip.text.x = element_text(face = "bold", size = 14, color = "grey20"),
    strip.background = element_rect(fill = "grey92", color = NA),
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
  file.path(tempdir(), "sc-celltype-sample-dodge-count-output")
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

# sc-celltype-grouped-horizontal-bar
#
# Clean, portable plotting example for Open Figure Modules. This linear R
# script reads only the synthetic CSV distributed with the module. External
# source images and unredistributed original code are not part of this module.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
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

# ---- 第三步：聚合成组水平比例 ----------------------------------------------
dat.plot <- counts %>%
  mutate(Sample = ifelse(as.character(Group) == "Treatment", "Treatment", "Control")) %>%
  group_by(Sample, celltype) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(value = 100 * n / sum(n)) %>%
  ungroup()
ctrl_order <- dat.plot %>%
  filter(Sample == "Control") %>%
  arrange(desc(value)) %>%
  pull(celltype)
dat.plot$Cell <- factor(dat.plot$celltype, levels = ctrl_order)
dat.plot$Sample <- factor(dat.plot$Sample, levels = c("Control", "Treatment"))
dat.plot$label <- paste0(round(dat.plot$value, 2), "%")

pal_fig1 <- c(
  "B-cells"             = "#E4D5EA",
  "cardiomyocyte cells" = "#6B3FA0",
  "NK-cells"            = "#C65A12",
  "Pericyte"            = "#9ECFC8",
  "T-cells"             = "#D4A017",
  "smooth muscle cells" = "#E39B1A",
  "myeloid cells"       = "#D46A8A",
  "Endothelial cells"   = "#3D8B4A",
  "Fibroblasts"         = "#6B7BB5"
)

# ---- 第四步：分面水平条 ----------------------------------------------------
p1 <- ggplot(dat.plot, aes(x = value, y = Cell)) +
  geom_col(aes(fill = Cell), width = 0.72) +
  geom_text(aes(x = value + 1.2, label = label), hjust = 0, size = 3, fontface = "bold") +
  facet_wrap(~ Sample) +
  scale_fill_manual(values = pal_fig1, guide = "none") +
  coord_cartesian(xlim = c(0, 45), clip = "off") +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "grey60", fill = "transparent"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10, face = "bold", colour = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 28, 8, 8)
  )

# ---- 最后一步：导出生成预览 ------------------------------------------------
output_args <- sub("^--output-dir=", "", args[grepl("^--output-dir=", args)])
output_env <- Sys.getenv("SFL_OUTPUT_DIR", unset = "")
out_dir <- if (length(output_args) && nzchar(output_args[[1]])) {
  output_args[[1]]
} else if (nzchar(output_env)) {
  output_env
} else {
  file.path(tempdir(), "sc-celltype-grouped-horizontal-bar-output")
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

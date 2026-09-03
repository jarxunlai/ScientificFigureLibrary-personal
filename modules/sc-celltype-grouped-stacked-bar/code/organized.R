# sc-celltype-grouped-stacked-bar
#
# Clean, portable plotting example for Open Figure Modules. This linear R
# script reads only the synthetic CSV distributed with the module. External
# source images and unredistributed original code are not part of this module.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})

WIDTH_IN <- 10
HEIGHT_IN <- 8

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

# ---- 第三步：组合计 --------------------------------------------------------
mydata <- counts %>%
  group_by(Group, celltype) %>%
  summarise(Quantity = sum(n), .groups = "drop")
ctrl_order <- mydata %>%
  filter(Group == "Control") %>%
  arrange(Quantity) %>%
  pull(celltype)
mydata$Celltype <- factor(mydata$celltype, levels = ctrl_order)
n_types <- nlevels(mydata$Celltype)
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_types)
names(my_colors) <- levels(mydata$Celltype)

# ---- 第四步：100% 堆叠柱 --------------------------------------------------
p1 <- ggplot(mydata, aes(x = Group, y = Quantity, fill = Celltype)) +
  geom_col(position = "fill", alpha = 0.7, width = 0.7) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Celltype", y = "Proportion", fill = "") +
  theme_bw(base_family = "serif") +
  theme(
    axis.title.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---- 最后一步：导出生成预览 ------------------------------------------------
output_args <- sub("^--output-dir=", "", args[grepl("^--output-dir=", args)])
output_env <- Sys.getenv("SFL_OUTPUT_DIR", unset = "")
out_dir <- if (length(output_args) && nzchar(output_args[[1]])) {
  output_args[[1]]
} else if (nzchar(output_env)) {
  output_env
} else {
  file.path(tempdir(), "sc-celltype-grouped-stacked-bar-output")
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

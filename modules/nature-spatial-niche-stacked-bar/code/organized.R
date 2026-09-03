# =============================================================================
# Nature 同款空转细胞生态位堆积柱
# =============================================================================
# 与 original.R：读合成比例表；ggplot2 4 用 position_fill(reverse=TRUE)
# 使 niche_9 在底、niche_1 在顶，对齐期刊原图。不是空转分析复现。
# =============================================================================

library(ggplot2)
library(readr)
library(dplyr)
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
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_font_family <- "Arial"

patient.prop <- read_csv(file.path(draft_dir, "data", "niche_proportions.csv"), show_col_types = FALSE)
sample_levels <- unique(patient.prop$patient)
niche_levels <- paste0("niche_", 1:9)
patient.prop$patient <- factor(patient.prop$patient, levels = sample_levels)
patient.prop$ctniche <- factor(patient.prop$ctniche, levels = rev(niche_levels))

cols <- c("#d51f26","#272e6a","#208a42","#89288f","#f47d2b","#fee500","#8a9fd1","#c06cab","#d8a767")
names(cols) <- niche_levels

p <- ggplot(patient.prop, aes(patient, Proportion, fill = ctniche)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.85) +
  xlab("") + ylab("cell type proportion") +
  scale_fill_manual(values = cols, breaks = niche_levels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic(base_size = 11, base_family = plot_font_family) +
  theme(
    axis.ticks.length = unit(0.15, "cm"),
    legend.position = "right",
    legend.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  guides(fill = guide_legend(title = NULL, ncol = 1, reverse = FALSE))

ggsave(file.path(out_dir, "niche_stacked_bar.png"), p, width = 10, height = 5.6, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "niche_stacked_bar.pdf"), p, width = 10, height = 5.6, device = cairo_pdf, bg = "white")
message("wrote niche_stacked_bar.png")

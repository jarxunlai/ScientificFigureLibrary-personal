# =============================================================================
# KEGG 变体条形图：蓝-黄-红，前 3 条白字
# =============================================================================
# 微信文 ggplot 第二版。数据与终图共用 kegg_top15.csv。
# =============================================================================

library(ggplot2)
library(readr)
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

dt <- as.data.frame(read_csv(file.path(draft_dir, "data", "kegg_top15.csv"), show_col_types = FALSE))
dt$Description <- factor(dt$Description, levels = rev(dt$Description))

p <- ggplot(data = dt, aes(x = -log10(p.adjust), y = Description, fill = -log10(p.adjust))) +
  scale_fill_gradientn(
    values = seq(0, 1, 0.1),
    colours = c("#4575b4", "#abd9e9", "#e0f3f8", "#ffffbf", "#fdae61", "#d73027", "#800026")
  ) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-Log10(pvalue)", y = "", title = "KEGG Pathway enrichment") +
  geom_text(
    size = 4.3, aes(x = 0.05, label = Description),
    hjust = 0, family = plot_font_family,
    color = rep(c("white", "black"), times = c(3, 12))
  ) +
  theme_classic(base_family = plot_font_family) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
  )

ggsave(file.path(out_dir, "kegg_bw_text_classic.png"), p, width = 7.2, height = 6.2, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "kegg_bw_text_classic.pdf"), p, width = 7.2, height = 6.2, device = cairo_pdf, bg = "white")
message("wrote kegg_bw_text_classic.png")

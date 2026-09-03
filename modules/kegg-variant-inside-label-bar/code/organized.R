# =============================================================================
# KEGG 变体条形图：通路名写在柱内
# =============================================================================
# 与 original.R：绑定 data/kegg_top15.csv；Seurat::NoLegend 改为 legend 位置。
# 15 条通路名对齐微信终图；p.adjust 为按终图色条合成，不是 pbmc3k 富集复现。
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
    colours = c("#ffe10f", "#fdae61", "#ee561a", "#d73027", "#dc0224"),
    name = "-log10(p.adjust)"
  ) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.9) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "", y = "", title = "") +
  geom_text(
    size = 4.3, aes(x = 0.05, label = Description),
    hjust = 0, family = plot_font_family,
    color = rep(c("white", "black"), times = c(7, 8))
  ) +
  theme_void(base_family = plot_font_family) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    legend.title = element_text(size = 13, angle = 90),
    legend.text = element_text(size = 11),
    legend.title.position = "left",
    legend.position = c(0.90, 0.15),
    plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
  )

ggsave(file.path(out_dir, "kegg_variant_bar.png"), p, width = 7.2, height = 6.2, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "kegg_variant_bar.pdf"), p, width = 7.2, height = 6.2, device = cairo_pdf, bg = "white")
message("wrote kegg_variant_bar.png")

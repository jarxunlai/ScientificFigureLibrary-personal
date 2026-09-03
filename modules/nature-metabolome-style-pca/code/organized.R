# =============================================================================
# Nature 风格四分组 PCA（椭圆 + 黑描边点）
# =============================================================================
# 作者 ggplot 美化段。pca_df 为合成坐标，方差比例用微信终图可见值
# PC1 33.16% / PC2 16.54%。不是 GSE126848 复现。
# 原文 p1 用 pca$x 叠一层点；这里 pca_df 已含 PC1/PC2，直接 shape=21。
# =============================================================================

library(ggplot2)
library(ggrepel)
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

pca_df <- as.data.frame(read_csv(file.path(draft_dir, "data", "pca_df.csv"), show_col_types = FALSE))
pca_df$group <- factor(pca_df$group, levels = c("healthy", "NAFLD", "NASH", "obese"))
var_explained <- c(0.3316, 0.1654)
mycol <- c(healthy = "#8ce5bb", NAFLD = "#7a7a7a", NASH = "#f6c6f4", obese = "#dcf7ea")

p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  stat_ellipse(alpha = 0.2, geom = "polygon", color = NA, type = "norm", level = 0.95) +
  geom_point(shape = 21, size = 3, stroke = 0.5, alpha = 0.9, color = "black") +
  geom_text_repel(aes(label = sample, color = group), size = 3, max.overlaps = 40, seed = 422, show.legend = FALSE) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  labs(
    x = paste0("PC1 (", round(var_explained[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 2), "%)")
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black", size = 12, family = plot_font_family),
    axis.text = element_text(color = "black", family = plot_font_family),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color = "black", size = 9, family = plot_font_family),
    legend.spacing.x = unit(0.06, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.background = element_blank(),
    legend.position = c(0.18, 0.82)
  )

ggsave(file.path(out_dir, "nature_style_pca.png"), p3, width = 6.4, height = 5.4, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "nature_style_pca.pdf"), p3, width = 6.4, height = 5.4, device = cairo_pdf, bg = "white")
message("wrote nature_style_pca.png")

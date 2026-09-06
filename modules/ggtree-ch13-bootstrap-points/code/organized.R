# =============================================================================
# Figure 13.2 分箱 bootstrap 圆点
# 数据：TDbook::text_RMI_tree
# =============================================================================

library(ggplot2)
library(ggtree)
library(treeio)
library(ragg)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nwk_text <- paste(readLines(file.path(root, "data", "RMI_tree.nwk"), warn = FALSE), collapse = "")
tree <- read.newick(text = nwk_text, node.label = "support")
root_id <- rootnode(tree)

p <- ggtree(tree, color = "black", size = 1.5, linetype = 1, right = TRUE) +
  geom_tiplab(size = 4.5, hjust = -0.060, fontface = "bold") +
  xlim(0, 0.09) +
  geom_point2(
    aes(subset = !isTip & node != root_id,
        fill = cut(support, c(0, 700, 900, 1000))),
    shape = 21, size = 4
  ) +
  theme_tree(legend.position = c(0.2, 0.2)) +
  scale_fill_manual(
    values = c("white", "grey", "black"),
    guide = "legend",
    name = "Bootstrap Percentage(BP)",
    breaks = c("(900,1e+03]", "(700,900]", "(0,700]"),
    labels = expression(BP >= 90, 70 <= BP * " < 90", BP < 70)
  )

ggsave(file.path(out_dir, "fig13_2_bootstrap_points.png"), p,
       width = 7.5, height = 8.6, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "fig13_2_bootstrap_points.pdf"), p,
       width = 7.5, height = 8.6, bg = "white", device = cairo_pdf)
file.copy(file.path(out_dir, "fig13_2_bootstrap_points.png"),
          file.path(root, "preview.png"), overwrite = TRUE)

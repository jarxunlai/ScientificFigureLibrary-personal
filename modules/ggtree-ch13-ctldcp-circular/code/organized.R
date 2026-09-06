# =============================================================================
# Figure 13.3 鸡 CTLDcp 环形树：Group II / V 高亮与禽特异扩张
# 数据：TDbook::tree_treenwk_30.4.19（Larsen et al. 2019）
# =============================================================================

library(ape)
library(ggplot2)
library(ggtree)
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

mytree <- read.tree(file.path(root, "data", "CTLDcp.nwk"))
tiplab <- mytree$tip.label
cls <- tiplab[grep("^ch", tiplab)]
labeltree <- groupOTU(mytree, cls)

p <- ggtree(labeltree, aes(color = group, linetype = group), layout = "circular") +
  scale_color_manual(values = c("#efad29", "#63bbd4")) +
  geom_nodepoint(color = "black", size = 0.1) +
  geom_tiplab(size = 2, color = "black")

p2 <- flip(p, 136, 110) %>%
  flip(141, 145) %>%
  ggtree::rotate(141) %>%
  ggtree::rotate(142) %>%
  ggtree::rotate(160) %>%
  ggtree::rotate(164) %>%
  ggtree::rotate(131)

dat <- data.frame(
  node = c(110, 88, 156, 136),
  fill = c("#229f8a", "#229f8a", "#229f8a", "#f9311f")
)
p3 <- p2 +
  geom_hilight(
    data = dat,
    mapping = aes(node = node, fill = I(fill)),
    alpha = 0.2,
    extendto = 1.4
  )

p4 <- p3 +
  geom_cladelab(
    node = 113,
    label = "Avian-specific expansion",
    align = TRUE, angle = -35, offset.text = 0.05,
    hjust = "center", fontsize = 2, offset = .2, barsize = .2
  )

p5 <- p4 +
  geom_nodelab(
    mapping = aes(
      x = branch,
      label = label,
      subset = !is.na(as.numeric(label)) & as.numeric(label) > 50
    ),
    size = 2, color = "black", nudge_y = 0.6
  )

p6 <- p5 +
  geom_cladelab(
    data = data.frame(node = c(114, 121), name = c("Subgroup A", "Subgroup B")),
    mapping = aes(node = node, label = name),
    align = TRUE, offset = .05, offset.text = .03, hjust = "center",
    barsize = .2, fontsize = 2, angle = "auto", horizontal = FALSE
  ) +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )

ggsave(file.path(out_dir, "fig13_3_ctldcp_circular.png"), p6,
       width = 7.5, height = 6.3, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "fig13_3_ctldcp_circular.pdf"), p6,
       width = 7.5, height = 6.3, bg = "white", device = cairo_pdf)
file.copy(file.path(out_dir, "fig13_3_ctldcp_circular.png"),
          file.path(root, "preview.png"), overwrite = TRUE)

# =============================================================================
# Figure 13.4 基因组位点结构 + Jaccard/BioNJ 树
# 数据：gggenes::example_genes（演示基因组）
# =============================================================================

library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)
library(ape)
library(readr)
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

example_genes <- read_csv(file.path(root, "data", "example_genes.csv"), show_col_types = FALSE)

get_genes <- function(data, genome) {
  filter(data, molecule == genome) %>% pull(gene)
}

g <- unique(example_genes[[1]])
n <- length(g)
d <- matrix(nrow = n, ncol = n)
rownames(d) <- colnames(d) <- g
genes <- lapply(g, get_genes, data = example_genes)

for (i in 1:n) {
  for (j in 1:i) {
    jaccard_sim <- length(intersect(genes[[i]], genes[[j]])) /
      length(union(genes[[i]], genes[[j]]))
    d[j, i] <- d[i, j] <- 1 - jaccard_sim
  }
}

tree <- ape::bionj(d)

# ggtree 4.0.5 的 geom_motif() 虽然有 label 参数，但内部仍 aes(label = label)。
# geom_facet 合并后 label 变成基因组名，箭头会被标成 Genome1 而不是 genA。
# 这里按书中意图：用 genE 对齐，文字用 gene 列。
geom_motif_gene <- function(mapping, data, on, align = "left", ...) {
  dd <- data[data$gene == on, ]
  mid <- dd$start + (dd$end - dd$start) / 2
  names(mid) <- dd$label
  adj <- unname(mid[data$label])
  data$start <- data$start - adj
  data$end <- data$end - adj
  ly_gene <- gggenes::geom_gene_arrow(
    mapping = aes(xmin = start, xmax = end, y = y, fill = gene),
    data = data, inherit.aes = FALSE, ...
  )
  ly_lab <- gggenes::geom_gene_label(
    mapping = aes(xmin = start, xmax = end, y = y, label = gene),
    data = data, align = align, inherit.aes = FALSE, ...
  )
  list(ly_gene, ly_lab)
}

p <- ggtree(tree, branch.length = "none") +
  geom_tiplab() +
  xlim_tree(5.5) +
  geom_facet(
    mapping = aes(xmin = start, xmax = end, fill = gene),
    data = example_genes, geom = geom_motif_gene, panel = "Alignment",
    on = "genE", align = "left"
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(strip.text = element_blank(), panel.spacing = unit(0, "cm"))

p <- facet_widths(p, widths = c(1, 2))

ggsave(file.path(out_dir, "fig13_4_genome_locus.png"), p,
       width = 9, height = 4, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "fig13_4_genome_locus.pdf"), p,
       width = 9, height = 4, bg = "white", device = cairo_pdf)
file.copy(file.path(out_dir, "fig13_4_genome_locus.png"),
          file.path(root, "preview.png"), overwrite = TRUE)

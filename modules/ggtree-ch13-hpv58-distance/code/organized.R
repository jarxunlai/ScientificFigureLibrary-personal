# =============================================================================
# Figure 13.1 HPV58 全基因组树 + 成对核苷酸距离
# organized：线性脚本 + 中文分节
# 数据：TDbook::tree_HPV58 / dna_HPV58_aln（Chen et al. 2017）
# 不重跑 GenBank 下载或 muscle 比对。
# =============================================================================

library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(tidyr)
library(tibble)
library(Biostrings)
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

# ---- 读入 HPV58 树与已比对的全基因组 --------------------------------------
tree <- read.tree(file.path(root, "data", "HPV58.nwk"))
clade <- c(A3 = 92, A1 = 94, A2 = 108, B1 = 156,
           B2 = 159, C = 163, D1 = 173, D2 = 176)
tree <- groupClade(tree, clade)
cols <- c(
  A1 = "#EC762F", A2 = "#CA6629", A3 = "#894418", B1 = "#0923FA",
  B2 = "#020D87", C = "#000000", D1 = "#9ACD32", D2 = "#08630A"
)

# ---- 树：谱系着色、clade bar、支持率 ----------------------------------------
p <- ggtree(tree, aes(color = group), ladderize = FALSE) %>%
  rotate(rootnode(tree)) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")),
              parse = TRUE, size = 2.5) +
  geom_treescale(x = 0, y = 1, width = 0.002) +
  scale_color_manual(
    values = c(cols, "black"),
    na.value = "black", name = "Lineage",
    breaks = c("A1", "A2", "A3", "B1", "B2", "C", "D1", "D2")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2(legend.position = c(.1, .88))

dat <- tibble(
  node = c(94, 108, 131, 92, 156, 159, 163, 173, 176, 172),
  name = c("A1", "A2", "A3", "A", "B1", "B2", "C", "D1", "D2", "D"),
  offset = c(0.003, 0.003, 0.003, 0.00315, 0.003, 0.003, 0.0031, 0.003, 0.003, 0.00315),
  offset.text = c(-.001, -.001, -.001, 0.0002, -.001, -.001, 0.0002, -.001, -.001, 0.0002),
  barsize = c(1.2, 1.2, 1.2, 2, 1.2, 1.2, 3.2, 1.2, 1.2, 2),
  extend = list(c(0, 0.5), 0.5, c(0.5, 0), 0, c(0, 0.5), c(0.5, 0), 0, c(0, 0.5), c(0.5, 0), 0)
) %>%
  dplyr::group_split(barsize)

p <- p +
  geom_cladelab(
    data = dat[[1]],
    mapping = aes(node = node, label = name, color = group,
                  offset = offset, offset.text = offset.text, extend = extend),
    barsize = 1.2, fontface = 3, align = TRUE
  ) +
  geom_cladelab(
    data = dat[[2]],
    mapping = aes(node = node, label = name,
                  offset = offset, offset.text = offset.text, extend = extend),
    barcolor = "darkgrey", textcolor = "darkgrey",
    barsize = 2, fontsize = 5, fontface = 3, align = TRUE
  ) +
  geom_cladelab(
    data = dat[[3]],
    mapping = aes(node = node, label = name,
                  offset = offset, offset.text = offset.text, extend = extend),
    barcolor = "darkgrey", textcolor = "darkgrey",
    barsize = 3.2, fontsize = 5, fontface = 3, align = TRUE
  ) +
  geom_strip(65, 71, "italic(B)", color = "darkgrey",
             offset = 0.00315, align = TRUE, offset.text = 0.0002,
             barsize = 2, fontsize = 5, parse = TRUE)

p <- p +
  geom_nodelab(aes(subset = (node == 92), label = "*"),
               color = "black", nudge_x = -.001, nudge_y = 1) +
  geom_nodelab(aes(subset = (node == 155), label = "*"),
               color = "black", nudge_x = -.0003, nudge_y = -1) +
  geom_nodelab(aes(subset = (node == 158), label = "95/92/1.00"),
               color = "black", nudge_x = -0.0001, nudge_y = -1, hjust = 1) +
  geom_nodelab(aes(subset = (node == 162), label = "98/97/1.00"),
               color = "black", nudge_x = -0.0001, nudge_y = -1, hjust = 1) +
  geom_nodelab(aes(subset = (node == 172), label = "*"),
               color = "black", nudge_x = -.0003, nudge_y = -1)

# ---- 成对 Hamming 距离（相对全长，百分数） --------------------------------
tl <- tree$tip.label
acc <- sub("\\w+\\|", "", tl)
names(tl) <- acc
tipseq_aln <- readDNAStringSet(file.path(root, "data", "HPV58_aln.fas"))
# Biostrings >= 2.77 把 stringDist 挪到 pwalign；这里用字符矩阵算 Hamming，不额外装包。
aln_mat <- as.matrix(tipseq_aln)
n_seq <- nrow(aln_mat)
n_pos <- ncol(aln_mat)
tipseq_d <- matrix(0, n_seq, n_seq, dimnames = list(names(tipseq_aln), names(tipseq_aln)))
for (i in seq_len(n_seq)) {
  for (j in seq_len(i)) {
    d_ij <- sum(aln_mat[i, ] != aln_mat[j, ]) / n_pos * 100
    tipseq_d[i, j] <- tipseq_d[j, i] <- d_ij
  }
}
dd <- as_tibble(tipseq_d)
dd$seq1 <- rownames(tipseq_d)
td <- gather(dd, seq2, dist, -seq1)
td$seq1 <- tl[td$seq1]
td$seq2 <- tl[td$seq2]
g <- p$data$group
names(g) <- p$data$label
td$clade <- g[td$seq2]

p2 <- p +
  geom_facet(
    panel = "Sequence Distance", data = td, geom = geom_point, alpha = .6,
    mapping = aes(x = dist, color = clade, shape = clade)
  ) +
  geom_facet(
    panel = "Sequence Distance", data = td, geom = geom_path, alpha = .6,
    mapping = aes(x = dist, group = seq2, color = clade)
  ) +
  scale_shape_manual(values = 1:8, guide = "none")

ggsave(file.path(out_dir, "fig13_1_hpv58_distance.png"), p2,
       width = 12, height = 12, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "fig13_1_hpv58_distance.pdf"), p2,
       width = 12, height = 12, bg = "white", device = cairo_pdf)
file.copy(file.path(out_dir, "fig13_1_hpv58_distance.png"),
          file.path(root, "preview.png"), overwrite = TRUE)

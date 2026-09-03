# =============================================================================
# ClusterGVis 风格：DEG 聚类热图 + 左侧趋势折线 + GO/KEGG 注释
# organized 版本：线性脚本 + 中文分节导航
# =============================================================================
# 这是绘图层复刻，不是原文分析复现。
# 微信文章未贴完整 R 代码；本脚本用 ComplexHeatmap 复刻 visCluster(
#   plotType = "both", lineSide = "left", markGenesSide = "left",
#   byGo/byKegg = "anno_block") 的版式。
# 表达矩阵为按参考图结构构造的合成 Z-score，不能声称复现 PRJNA975358
# 或 DESeq2/clusterProfiler 结果。
# =============================================================================

library(dplyr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)
library(ragg)

# ---- 路径：以本脚本所在目录为基准 ------------------------------------------
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      normalizePath(getwd())
    }
  }
)
if (basename(script_dir) == "code") {
  draft_dir <- dirname(script_dir)
} else {
  draft_dir <- script_dir
}

data_dir <- file.path(draft_dir, "data")
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_png <- file.path(out_dir, "clustergvis_deg_go_kegg.png")
output_pdf <- file.path(out_dir, "clustergvis_deg_go_kegg.pdf")
preview_png <- file.path(draft_dir, "preview.png")

plot_font_family <- "Arial"
figure_width_in <- 12.8
figure_height_in <- 9.2
figure_dpi <- 300

# ---- 第一步：固定样本、聚类顺序和标记基因 ----------------------------------
# 文章思考过程写 vst_counts.tsv 有 14 个样本。参考图顶部为 WT_control /
# KRAS_G12D 两组；第三对样本名在截图里不够清晰，这里按 E262 记录，
# 并保留 E95_wt / E95_g12d，使列数为 14。
# 行模块按参考图从上到下 C3 → C1 → C4 → C2，对应 visCluster 的
# clusterOrder = c(3, 1, 4, 2)。
wt_samples <- c(
  "E214_wt", "E218_wt", "E262_wt", "E263_wt", "E266_wt", "E9_wt", "E95_wt"
)
kras_samples <- c(
  "E214_g12d", "E218_g12d", "E262_g12d", "E263_g12d",
  "E266_g12d", "E9_g12d", "E95_g12d"
)
sample_ids <- c(wt_samples, kras_samples)
sample_group <- factor(
  c(rep("WT_control", length(wt_samples)), rep("KRAS_G12D", length(kras_samples))),
  levels = c("WT_control", "KRAS_G12D")
)
names(sample_group) <- sample_ids

cluster_sizes <- c(C3 = 41L, C1 = 50L, C4 = 39L, C2 = 150L)
cluster_order <- names(cluster_sizes)

mark_genes <- c(
  C3 = list(c("Col12a1", "Aspn", "Fbn2")),
  C1 = list(character(0)),
  C4 = list(c("H2-D1", "B2m", "Tap2")),
  C2 = list(c(
    "Cxcl9", "Oas3", "H2-K1", "Ccl5", "H2-Q7", "Ifit2",
    "Cxcl10", "Zbp1", "Tap1", "Oas1b", "Mx2", "Oas2",
    "Ifit1", "Isg15", "Ddx58", "Stat2", "Oas1a", "Irf7",
    "Rsad2", "Ifih1", "Stat1", "Mx1", "Tlr3", "Tapbp"
  ))
)

# 参考图：C3 绿、C1 青、C4 浅蓝、C2 靛蓝。
ct_anno_col <- c(
  C3 = "#7CB342",
  C1 = "#26A69A",
  C4 = "#4FC3F7",
  C2 = "#1A237E"
)
sample_col <- c(
  WT_control = "#6B87B3",
  KRAS_G12D = "#C27B8A"
)
ht_col <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#08519C", "white", "#A50F15")
)

# ---- 第二步：构造与参考图趋势一致的合成 Z-score -----------------------------
# 每个 cluster 给一条 14 点中位轮廓，再叠加基因噪声。
# C3/C1：WT 高、KRAS 低；C4/C2：WT 低、KRAS 高；C2 最后一个 KRAS 样本回落。
set.seed(20260701)

cluster_profile <- list(
  C3 = c(1.55, 1.35, 1.25, 1.05, 0.70, 0.25, 0.05, -0.85, -1.05, -1.15, -1.25, -1.10, -1.20, -1.05),
  C1 = c(1.45, 1.10, 1.30, 0.85, 1.00, 0.40, 0.15, -0.55, -0.95, -0.70, -1.05, -1.15, -0.80, -1.00),
  C4 = c(-1.15, -1.05, -0.95, -0.85, -0.75, -0.60, -0.35, 0.55, 1.05, 1.35, 1.45, 1.20, 1.05, 0.85),
  C2 = c(-1.35, -1.15, -1.25, -1.05, -1.20, -0.95, -0.75, 0.75, 1.25, 1.55, 1.40, 1.30, 1.05, 0.15)
)

make_cluster_genes <- function(cluster_id, n, marked) {
  n_marked <- length(marked)
  extra <- if (n > n_marked) {
    sprintf("%s_g%02d", cluster_id, seq_len(n - n_marked))
  } else {
    character(0)
  }
  genes <- c(marked, extra)
  if (length(genes) > n) {
    genes <- genes[seq_len(n)]
  }
  profile <- cluster_profile[[cluster_id]]
  mat <- vapply(seq_along(genes), function(i) {
    jitter_sd <- if (genes[i] %in% marked) 0.18 else 0.32
    val <- profile + rnorm(length(profile), mean = 0, sd = jitter_sd)
    pmin(pmax(val, -2.2), 2.2)
  }, numeric(length(profile)))
  mat <- t(mat)
  colnames(mat) <- sample_ids
  rownames(mat) <- genes
  as.data.frame(mat, check.names = FALSE) |>
    tibble::rownames_to_column("gene") |>
    mutate(cluster = cluster_id, .before = 2)
}

wide_list <- lapply(cluster_order, function(cid) {
  make_cluster_genes(cid, cluster_sizes[[cid]], mark_genes[[cid]])
})
wide_res <- bind_rows(wide_list)

# 长词条在词边界手工换行，避免 textbox_grob 的 word_wrap 把空格吃掉。
go_terms <- tibble::tribble(
  ~id, ~term,
  "C3", "positive regulation of canonical\nWnt signaling pathway",
  "C3", "osteoblast differentiation",
  "C3", "regulation of epithelial\ncell proliferation",
  "C3", "regulation of transmembrane receptor\nprotein serine/threonine kinase\nsignaling pathway",
  "C3", "positive regulation of Wnt\nsignaling pathway",
  "C1", "extracellular matrix organization",
  "C1", "homophilic cell adhesion via\nplasma membrane adhesion molecules",
  "C1", "external encapsulating\nstructure organization",
  "C1", "extracellular structure organization",
  "C1", "cell-cell adhesion via\nplasma-membrane adhesion molecules",
  "C4", "regulation of leukocyte\nmediated cytotoxicity",
  "C4", "regulation of cell killing",
  "C4", "antigen processing and presentation\nof peptide antigen via MHC class Ib",
  "C4", "cell killing",
  "C4", "positive regulation of leukocyte\nmediated cytotoxicity",
  "C2", "response to virus",
  "C2", "defense response to virus",
  "C2", "response to interferon-beta",
  "C2", "cellular response to interferon-beta",
  "C2", "activation of innate immune response"
)

kegg_terms <- tibble::tribble(
  ~id, ~term,
  "C3", "Axon guidance",
  "C3", "Calcium signaling pathway",
  "C3", "Breast cancer",
  "C3", "Wnt signaling pathway",
  "C3", "Apelin signaling pathway",
  "C1", "Protein digestion and absorption",
  "C1", "Hippo signaling pathway -\nmultiple species",
  "C1", "Cytoskeleton in muscle cells",
  "C1", "Cadherin signaling",
  "C4", "Antigen processing and presentation",
  "C4", "Herpes simplex virus 1 infection",
  "C4", "Human cytomegalovirus infection",
  "C4", "IL-17 signaling pathway",
  "C4", "Epstein-Barr virus infection",
  "C2", "Herpes simplex virus 1 infection",
  "C2", "Influenza A",
  "C2", "Hepatitis C",
  "C2", "NOD-like receptor signaling pathway",
  "C2", "Epstein-Barr virus infection"
)

# 每条 term 固定配色，避免 circlize::rand_color 每次重跑都变。
go_col <- c(
  "#43A047", "#00897B", "#1E88E5", "#8E24AA", "#D81B60",
  "#7CB342", "#00ACC1", "#5E35B1", "#3949AB", "#F06292",
  "#66BB6A", "#26C6DA", "#42A5F5", "#7E57C2", "#FFA726",
  "#2E7D32", "#FB8C00", "#8E24AA", "#EC407A", "#6D4C41"
)
kegg_col <- c(
  "#7CB342", "#26A69A", "#E53935", "#FB8C00", "#F9A825",
  "#8E24AA", "#43A047", "#5C6BC0", "#EF6C00",
  "#26A69A", "#7E57C2", "#EC407A", "#FFA726", "#8D6E63",
  "#AB47BC", "#42A5F5", "#EC407A", "#7E57C2", "#8D6E63"
)

# ---- 第三步：写出可检查的合成表 --------------------------------------------
write_csv(wide_res, file.path(data_dir, "synthetic_zscore.csv"))
write_csv(go_terms, file.path(data_dir, "go_terms.csv"))
write_csv(kegg_terms, file.path(data_dir, "kegg_terms.csv"))

# ---- 第四步：整理 visCluster 同款矩阵与行切分 -------------------------------
mat <- as.matrix(wide_res[, sample_ids])
rownames(mat) <- wide_res$gene
subgroup <- factor(wide_res$cluster, levels = cluster_order)
align_to <- split(seq_len(nrow(mat)), subgroup)

# ---- 第五步：顶部样本分组、右侧 cluster 色条、左侧标记基因 ------------------
topanno <- HeatmapAnnotation(
  sample = sample_group,
  col = list(sample = sample_col),
  gp = gpar(col = "white"),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    sample = list(title = "sample", nrow = 2)
  )
)

cluster_id_anno <- anno_block(
  align_to = align_to,
  panel_fun = function(index, nm) {
    grid.text(
      nm,
      rot = 90,
      gp = gpar(
        fontsize = 10,
        fontface = "bold",
        fontfamily = plot_font_family
      )
    )
  },
  width = unit(4.5, "mm"),
  which = "row"
)

anno_block <- anno_block(
  align_to = align_to,
  panel_fun = function(index, nm) {
    fill_col <- ct_anno_col[[nm]]
    txt_col <- if (nm == "C2") "white" else "grey15"
    grid.rect(gp = gpar(fill = fill_col, col = NA))
    grid.text(
      label = paste("Num:", length(index)),
      rot = 90,
      gp = gpar(
        col = txt_col,
        fontsize = 8,
        fontfamily = plot_font_family
      )
    )
  },
  width = unit(4.2, "mm"),
  which = "row"
)

all_mark <- unlist(mark_genes[cluster_order], use.names = FALSE)
mark_col <- setNames(
  ct_anno_col[wide_res$cluster[match(all_mark, wide_res$gene)]],
  all_mark
)
gene_mark <- anno_mark(
  at = match(all_mark, rownames(mat)),
  labels = all_mark,
  which = "row",
  side = "left",
  labels_gp = gpar(
    fontface = "italic",
    fontsize = 6.6,
    col = unname(mark_col),
    fontfamily = plot_font_family
  ),
  link_width = unit(2.8, "mm"),
  padding = unit(0.4, "mm")
)

# ---- 第六步：左侧 cluster 中位趋势（visCluster panel_fun 的线性写法） --------
# 中位线颜色与 visCluster 默认 mlineCol = "#CC3333" 一致。
# y 轴用全局矩阵范围，和 visCluster 一样，这样各 cluster 的高低可比较。
rg <- range(mat)
line_anno <- anno_link(
  align_to = align_to,
  which = "row",
  side = "left",
  size = unit(1.35, "cm"),
  gap = unit(1.6, "mm"),
  width = unit(2.05, "cm"),
  link_gp = gpar(fill = "grey90", col = NA),
  panel_fun = function(index, nm) {
    pushViewport(viewport(xscale = c(-0.1, 1.1), yscale = c(0, 1)))
    grid.rect(gp = gpar(fill = "white", col = "grey35", lwd = 0.8))
    tmpmat <- mat[index, , drop = FALSE]
    mdia <- apply(tmpmat, 2, median)
    grid.lines(
      x = rescale(seq_len(ncol(tmpmat)), to = c(0.08, 0.92)),
      y = rescale(mdia, to = c(0.08, 0.72), from = c(rg[1] - 0.5, rg[2] + 0.5)),
      default.units = "native",
      gp = gpar(lwd = 2.0, col = "#CC3333")
    )
    grid.rect(
      x = 0.50,
      y = 0.88,
      width = 0.96,
      height = 0.20,
      default.units = "native",
      gp = gpar(fill = "white", col = "grey40", lwd = 0.55)
    )
    grid.text(
      paste("Gene size:", length(index)),
      x = 0.50,
      y = 0.88,
      default.units = "native",
      gp = gpar(
        fontsize = 5.6,
        fontface = "italic",
        col = "black",
        fontfamily = plot_font_family
      )
    )
    popViewport()
  }
)

# ---- 第七步：GO / KEGG 文本盒（anno_textbox + anno_block） -------------------
# visCluster 把 enrichCluster 表的 id/term 列喂给 anno_textbox。
# C1 的 KEGG 只有 4 条，不能按 cluster_num * 5 生成颜色。
make_term_list <- function(term_df, colors) {
  term_df <- term_df |>
    mutate(col = colors, fontsize = 6.2, fontfamily = plot_font_family)
  out <- lapply(cluster_order, function(cid) {
    tmp <- term_df[term_df$id == cid, , drop = FALSE]
    data.frame(
      text = tmp$term,
      col = tmp$col,
      fontsize = tmp$fontsize,
      fontfamily = tmp$fontfamily,
      stringsAsFactors = FALSE
    )
  })
  names(out) <- cluster_order
  out
}

go_list <- make_term_list(go_terms, go_col)
kegg_list <- make_term_list(kegg_terms, kegg_col)

go_box <- anno_textbox(
  align_to,
  go_list,
  word_wrap = FALSE,
  add_new_line = TRUE,
  side = "right",
  background_gp = gpar(fill = "grey95", col = "grey50"),
  by = "anno_block",
  max_width = unit(62, "mm"),
  padding = unit(3, "pt"),
  line_space = unit(2, "pt")
)
kegg_box <- anno_textbox(
  align_to,
  kegg_list,
  word_wrap = FALSE,
  add_new_line = TRUE,
  side = "right",
  background_gp = gpar(fill = "grey95", col = "grey50"),
  by = "anno_block",
  max_width = unit(48, "mm"),
  padding = unit(3, "pt"),
  line_space = unit(2, "pt")
)

# ---- 第八步：拼 Heatmap 并导出 ----------------------------------------------
ht_opt(message = FALSE)

ht <- Heatmap(
  matrix = mat,
  name = "Z-score",
  col = ht_col,
  border = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  column_split = sample_group,
  column_gap = unit(2.2, "mm"),
  row_split = subgroup,
  row_gap = unit(1.4, "mm"),
  row_title = NULL,
  column_title_gp = gpar(fontsize = 10, fontfamily = plot_font_family),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 7.5, fontfamily = plot_font_family),
  column_names_side = "top",
  top_annotation = topanno,
  left_annotation = rowAnnotation(
    cluster_id = cluster_id_anno,
    line = line_anno,
    gene = gene_mark,
    annotation_name_gp = gpar(fontsize = 0),
    gap = unit(1.2, "mm")
  ),
  right_annotation = rowAnnotation(
    cluster = anno_block,
    GO = go_box,
    KEGG = kegg_box,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 8, fontfamily = plot_font_family),
    gap = unit(1.0, "mm")
  ),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(2.4, "cm"),
    grid_width = unit(3.5, "mm")
  ),
  width = unit(7.6, "cm"),
  use_raster = FALSE
)

draw_figure <- function() {
  draw(
    ht,
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(1.5, 1.5, 1.5, 1.5), "mm")
  )
}

agg_png(
  output_png,
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = figure_dpi,
  background = "white"
)
draw_figure()
dev.off()

cairo_pdf(
  output_pdf,
  width = figure_width_in,
  height = figure_height_in,
  onefile = FALSE
)
draw_figure()
dev.off()

file.copy(output_png, preview_png, overwrite = TRUE)

message("wrote ", output_png)
message("wrote ", output_pdf)

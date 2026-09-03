# ---- 生成高保真合成示例数据（仅供验证绘图代码可执行性） ----------------------
# 目标：让 original.R 跑出的图在结构、密度、标签上尽量贴近论文 panel C 截图。
# 注意：模块基因数必须与 original.R 中硬编码标注一致（up=438, down=416），
#       GO term 直接采用截图中可见的 10 条条目原文。

set.seed(42)

out_dir <- "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 细胞类型：名称与分组取自截图底部轴标签 --------------------------------
cell_anno <- tibble::tribble(
  ~tissue,         ~cell_type,
  "Hypothalamus",  "MGC",
  "Hypothalamus",  "MOC",
  "Hypothalamus",  "Astro",
  "Hypothalamus",  "Pericytes",
  "Hypothalamus",  "MAC",
  "Hypothalamus",  "T cells",
  "Pituitary",     "Gona",
  "Pituitary",     "Cort",
  "Pituitary",     "Lac",
  "Pituitary",     "Som",
  "Pituitary",     "PSC_Fol",
  "Pituitary",     "Tany",
  "Pituitary",     "TSC",
  "Pituitary",     "Endo",
  "Pituitary",     "MAC.2",
  "Pituitary",     "T cells.2",
  "Ovary",         "GC",
  "Ovary",         "SC",
  "Ovary",         "Theca",
  "Ovary",         "LUM+ SC",
  "Ovary",         "LUM- SC",
  "Ovary",         "SMC",
  "Ovary",         "Pericytes.2",
  "Ovary",         "Endo.2",
  "Ovary",         "FBC",
  "Ovary",         "MAC.3",
  "Ovary",         "T cells.3",
  "Uterus",        "Myo1",
  "Uterus",        "Myo2",
  "Uterus",        "Pericytes.3",
  "Uterus",        "SFRP4+ SC",
  "Uterus",        "LUM+ SC.2",
  "Uterus",        "Endo.3",
  "Uterus",        "FBC.2",
  "Uterus",        "MAC.4",
  "Uterus",        "T cells.4"
) |>
  dplyr::mutate(cell_index = dplyr::row_number())

readr::write_csv(cell_anno, file.path(out_dir, "input_cell_annotation.csv"))

# ---- 基因模块：行数匹配 original.R 的硬编码标注 438 / 416 ---------------------
# 参考图底部还有一段未标模块名的稀疏区域；用 300 个 Other 基因模拟，使高度比例约为
# Common up : Common down : Other = 38% : 36% : 26%。
n_up <- 438L
n_down <- 416L
n_other <- 300L

# 固定宽度编号保证按字符排序后仍是 001, 002, ...；否则 GeneUp1/GeneUp10 的
# 字典序会打乱我们设计的纵向密度梯度。
genes <- tibble::tibble(
  gene = c(
    sprintf("GeneUp%03d", seq_len(n_up)),
    sprintf("GeneDown%03d", seq_len(n_down)),
    sprintf("GeneOther%03d", seq_len(n_other))
  ),
  module = rep(c("Common up", "Common down", "Other"), c(n_up, n_down, n_other))
)

# ---- 差异状态矩阵 ------------------------------------------------------------
# 参考图的宏观纹理：
#   - Common up：上部粉色较密、向下逐渐稀疏；Ovary 最密，Uterus 最稀。
#   - Common down：上部至中部蓝色较密、向下衰减；Uterus 最密，Hypothalamus 最稀。
#   - Other：大面积白色，仅有零星粉/蓝短条。
# 同时引入行趋势、组织趋势、细胞列趋势和少量随机噪声，避免每列密度完全相同。

sigmoid <- function(x) 1 / (1 + exp(-x))

build_module_rows <- function(gene_vec, module_name, n_in_module) {
  idx_map <- setNames(seq_along(gene_vec), gene_vec)
  match_dir <- if (module_name == "Common up") "Upregulated in T1DM" else "Downregulated in T1DM"
  opp_dir <- if (module_name == "Common up") "Downregulated in T1DM" else "Upregulated in T1DM"

  grid <- expand.grid(
    gene = gene_vec, cell_type = cell_anno$cell_type,
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::left_join(genes, by = "gene") |>
    dplyr::left_join(cell_anno, by = "cell_type")

  # 0 = 模块顶部，1 = 模块底部（与 original.R 的 rev(gene_order) 一致）。
  row_frac <- (idx_map[grid$gene] - 1) / (n_in_module - 1)
  col_frac <- (grid$cell_index - 1) / (max(cell_anno$cell_index) - 1)

  if (module_name == "Common up") {
    tissue_mult <- c(Hypothalamus = 0.78, Pituitary = 1.05, Ovary = 1.18, Uterus = 0.68)
    col_wave <- 0.80 + 0.22 * sin(2 * pi * (col_frac * 2.2 + 0.10))^2
    p_match <- (0.04 + 0.78 * sigmoid(8.5 * (0.53 - row_frac))) *
      unname(tissue_mult[grid$tissue]) * col_wave
    p_opp <- 0.003 + 0.004 * row_frac
  } else {
    tissue_mult <- c(Hypothalamus = 0.62, Pituitary = 0.92, Ovary = 1.00, Uterus = 1.18)
    col_wave <- 0.82 + 0.24 * cos(2 * pi * (col_frac * 2.0 + 0.18))^2
    p_match <- (0.025 + 0.88 * sigmoid(9.0 * (0.58 - row_frac))) *
      unname(tissue_mult[grid$tissue]) * col_wave
    p_opp <- 0.003 + 0.010 * sigmoid(8 * (row_frac - 0.82))
  }

  p_match <- pmin(p_match, 0.96)
  p_opp <- pmin(p_opp, 0.025)
  u <- runif(nrow(grid))
  grid$regulation <- ifelse(
    u < p_match, match_dir,
    ifelse(u < p_match + p_opp, opp_dir, "Unchanged")
  )

  grid |>
    dplyr::select(gene, module, cell_type, tissue, regulation)
}

deg_up <- build_module_rows(
  genes$gene[genes$module == "Common up"], "Common up", n_up
)
deg_down <- build_module_rows(
  genes$gene[genes$module == "Common down"], "Common down", n_down
)

# Other 模块：用连续行段制造参考图中稀疏的水平短条，而不是完全独立的像素噪声。
deg_other <- expand.grid(
  gene = genes$gene[genes$module == "Other"],
  cell_type = cell_anno$cell_type,
  stringsAsFactors = FALSE
) |>
  tibble::as_tibble() |>
  dplyr::mutate(module = "Other") |>
  dplyr::left_join(cell_anno, by = "cell_type") |>
  dplyr::mutate(regulation = "Unchanged")

set.seed(43)
other_genes <- unique(deg_other$gene)
for (i in seq_len(34)) {
  g <- sample(other_genes, 1)
  start <- sample(seq_len(nrow(cell_anno) - 3), 1)
  len <- sample(2:8, 1)
  direction <- sample(
    c("Upregulated in T1DM", "Downregulated in T1DM"), 1,
    prob = c(0.58, 0.42)
  )
  hit_cells <- cell_anno$cell_type[start:min(start + len - 1, nrow(cell_anno))]
  deg_other$regulation[deg_other$gene == g & deg_other$cell_type %in% hit_cells] <- direction
}

deg <- dplyr::bind_rows(
  deg_up,
  deg_down,
  deg_other |> dplyr::select(gene, module, cell_type, tissue, regulation)
)
readr::write_csv(deg, file.path(out_dir, "input_deg_heatmap.csv"))

# ---- GO 富集：条目与柱长按截图读取 ------------------------------------------
go <- tibble::tribble(
  ~module,        ~term,                                          ~neg_log10_p,
  "Common up",    "oxidative phosphorylation",                    23.4,
  "Common up",    "positive regulation of leukocyte activation",   4.2,
  "Common up",    "response to steroid hormone",                   4.0,
  "Common up",    "positive regulation of apoptotic process",      2.5,
  "Common up",    "reproductive system development",               2.2,
  "Common down",  "blood vessel development",                     18.3,
  "Common down",  "negative regulation of cell migration",        14.8,
  "Common down",  "regulation of proteolysis",                    14.1,
  "Common down",  "response to endoplasmic reticulum stress",     14.0,
  "Common down",  "regulation of endothelial cell proliferation",  5.7
)

readr::write_csv(go, file.path(out_dir, "input_go_enrichment.csv"))

cat("genes per module:", n_up, n_down, n_other, "\n")
cat("cell types:", nrow(cell_anno), "\n")
cat("deg rows:", nrow(deg), "\n")

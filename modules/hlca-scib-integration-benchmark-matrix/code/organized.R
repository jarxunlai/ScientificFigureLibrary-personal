# =============================================================================
# HLCA Supplementary Figure 1：scIB 整合基准总表
# organized 版本：线性脚本 + 中文分节导航
# =============================================================================
# 与 source/official_scripts 的关系：
#   计分、排序、几何和 ColorBrewer 色阶跟 Lisa Sikkema 的
#   plotSingleTaskRNA.R + knit_table.R 一致（batch 权重 0.4）。
# 与官方脚本的可见差异：
#   1. 路径相对本草稿，支持 Rscript / source() / 仓库根目录，不再写死 macOS 路径；
#   2. 不依赖 dynutils / Hmisc / plyr / ggimage / cowplot / png；
#   3. Output 图标改为 ggplot 矢量，不读官方 PNG；
#   4. ggplot2 4 用 linewidth；
#   5. 补上期刊图里的星号脚注（官方 knit_table 没有这段文字）；
#   6. 补上期刊补图页左上角 “SUPPLEMENTARY FIGURES” 标题。
# 数据：HLCA_reproducibility 仓库的 metrics_scgen_added.csv（scIB 0.1.1 输出）。
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(stringr)
library(scales)
library(grid)
library(ragg)

# ---- 路径：以本脚本所在目录为基准 ------------------------------------------
# 需要能找到 data/metrics_scgen_added.csv。下面几种启动方式都应工作：
#   Rscript code/organized.R
#   source("drafts/.../code/organized.R")  # 从仓库根
#   工作目录就在本草稿内
candidates <- character()
ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!is.null(ofile) && nzchar(ofile)) {
  candidates <- c(candidates, dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
}
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) && nzchar(file_arg[[1]])) {
  candidates <- c(
    candidates,
    dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = FALSE))
  )
}
candidates <- c(
  candidates,
  file.path(getwd(), "code"),
  file.path(getwd(), "drafts", "hlca-scib-integration-benchmark-matrix", "code"),
  getwd()
)
draft_dir <- NULL
for (cand in unique(candidates)) {
  parent <- if (identical(basename(cand), "code")) dirname(cand) else cand
  if (file.exists(file.path(parent, "data", "metrics_scgen_added.csv"))) {
    draft_dir <- normalizePath(parent, winslash = "/", mustWork = TRUE)
    break
  }
}
if (is.null(draft_dir)) {
  stop("找不到 data/metrics_scgen_added.csv；请从仓库根目录或本草稿目录运行。")
}

metrics_csv <- file.path(draft_dir, "data", "metrics_scgen_added.csv")
official_summary_csv <- file.path(draft_dir, "data", "optional_reference_summary.csv")
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_font_family <- "Arial"
weight_batch <- 0.4

# ---- 第一步：读 scIB 指标宽表 ----------------------------------------------
# 每行 = 一种方法 × 一种预处理 × 一种输出。空字段是该输出类型算不了的指标。
metrics_raw <- read.csv(metrics_csv, check.names = FALSE, stringsAsFactors = FALSE)
names(metrics_raw)[1] <- "run_id"

# ---- 第二步：把列名改成期刊图上的指标名 ------------------------------------
metrics <- names(metrics_raw)[-1]
metrics <- chartr(".", intToUtf8(47), metrics)
metrics <- gsub("_", " ", metrics)
rename_from <- c(
  "ASW label", "ASW label/batch", "cell cycle conservation", "hvg overlap",
  "trajectory", "graph conn", "iLISI", "cLISI"
)
rename_to <- c(
  "Cell type ASW", "Batch ASW", "CC conservation", "HVG conservation",
  "trajectory conservation", "graph connectivity", "graph iLISI", "graph cLISI"
)
hit <- match(metrics, rename_from)
metrics[!is.na(hit)] <- rename_to[hit[!is.na(hit)]]

group_batch <- c("PCR batch", "Batch ASW", "graph iLISI", "graph connectivity", "kBET")
group_bio <- c(
  "NMI cluster/label", "ARI cluster/label", "Cell type ASW",
  "isolated label F1", "isolated label silhouette", "graph cLISI",
  "CC conservation", "HVG conservation", "trajectory conservation"
)
n_metrics_batch_original <- sum(group_batch %in% metrics)
n_metrics_bio_original <- sum(group_bio %in% metrics)

matching_order <- match(c(group_batch, group_bio), metrics)
metrics_ord <- metrics[matching_order[!is.na(matching_order)]]

# ---- 第三步：从 run_id 拆方法 / 预处理 / 输出 ------------------------------
# 官方路径形如 /lung_atlas_fixed/metrics/{scaled|unscaled}/{hvg|full_feature}/{method}_{output}
methods_info_full <- as.character(metrics_raw$run_id)
methods_info_full <- sub("^/", "", methods_info_full)

parts <- str_split(methods_info_full, "/", simplify = TRUE)
scaling <- parts[, 3]
hvg <- ifelse(parts[, 4] == "hvg", "HVG", "FULL")
method_token <- parts[, 5]
methods_key <- sub("_(knn|embed|full)$", "", method_token)
method_groups <- str_extract(method_token, "(knn|embed|full)$")
method_groups <- recode(
  method_groups,
  knn = "graph",
  embed = "embed",
  full = "gene",
  .default = method_groups
)

method_label <- recode(
  methods_key,
  seurat = "Seurat v3 CCA",
  seuratrpca = "Seurat v3 RPCA",
  mnn = "MNN",
  bbknn = "BBKNN",
  trvae = "trVAE",
  scvi = "scVI",
  liger = "LIGER",
  combat = "ComBat",
  saucie = "SAUCIE",
  fastmnn = "fastMNN",
  desc = "DESC",
  scanvi = "scANVI*",
  scgen = "scGen*",
  scanorama = "Scanorama",
  harmony = "Harmony",
  conos = "Conos",
  .default = methods_key
)

# ---- 第四步：丢掉全空列，按 scIB 规则算三类总分 ----------------------------
# 个体指标先按列 min-max 到 [0,1]，再对 batch / bio 组内取均值。
# Overall = 0.4 * batch + 0.6 * bio。图上圆点仍用未缩放的原始分数。
metrics_tab <- metrics_raw[, -1, drop = FALSE]
metrics_tab[metrics_tab == ""] <- NA
names(metrics_tab) <- metrics
metrics_tab <- cbind(
  Method = method_label,
  metrics_tab[, metrics_ord, drop = FALSE],
  stringsAsFactors = FALSE
)

na_col <- vapply(
  metrics_tab,
  function(x) sum(is.na(x)) == nrow(metrics_tab),
  logical(1)
)
n_metrics_batch <- n_metrics_batch_original - sum(names(metrics_tab)[na_col] %in% group_batch)
n_metrics_bio <- n_metrics_bio_original - sum(names(metrics_tab)[na_col] %in% group_bio)
metrics_tab <- metrics_tab[, !na_col, drop = FALSE]

scale_minmax <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || rng[1] == rng[2]) {
    return(ifelse(is.na(x), NA_real_, 0.5))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

scaled_metrics <- as.matrix(metrics_tab[, -1, drop = FALSE])
scaled_metrics <- apply(scaled_metrics, 2, scale_minmax)
score_group_batch <- rowMeans(scaled_metrics[, seq_len(n_metrics_batch), drop = FALSE], na.rm = TRUE)
score_group_bio <- rowMeans(
  scaled_metrics[, (n_metrics_batch + 1):ncol(scaled_metrics), drop = FALSE],
  na.rm = TRUE
)
score_all <- weight_batch * score_group_batch + (1 - weight_batch) * score_group_bio

metrics_tab <- cbind(
  Method = metrics_tab$Method,
  Output = method_groups,
  Features = hvg,
  Scaling = scaling,
  `Overall Score` = score_all,
  `Batch Correction` = score_group_batch,
  metrics_tab[, seq_len(n_metrics_batch) + 1, drop = FALSE],
  `Bio conservation` = score_group_bio,
  metrics_tab[, (n_metrics_batch + 2):ncol(metrics_tab), drop = FALSE],
  stringsAsFactors = FALSE
)

metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`, decreasing = TRUE), ]
rownames(metrics_tab) <- NULL

summary_out <- file.path(out_dir, "lung_atlas_fixed_summary_scores_recomputed.csv")
write.csv(metrics_tab, summary_out, row.names = FALSE, quote = FALSE)

rows_na <- which(is.na(metrics_tab$`Overall Score`))
if (length(rows_na)) {
  metrics_tab <- metrics_tab[-rows_na, ]
}

# ---- 第五步：和官方 summary 核对总分，防止改了计分还继续画图 --------------
if (file.exists(official_summary_csv)) {
  official <- read.csv(official_summary_csv, check.names = FALSE, stringsAsFactors = FALSE)
  official <- official[!is.na(official$`Overall Score`), ]
  stopifnot(nrow(official) == nrow(metrics_tab))
  stopifnot(identical(official$Method, metrics_tab$Method))
  stopifnot(identical(official$Output, metrics_tab$Output))
  stopifnot(identical(official$Features, metrics_tab$Features))
  stopifnot(max(abs(official$`Overall Score` - metrics_tab$`Overall Score`)) < 1e-8)
}

# ---- 第六步：行、列几何（官方 knit_table 的间距） --------------------------
row_height <- 1.1
row_space <- 0.1
col_width <- 1.1
col_space <- 0.2
col_bigspace <- 0.5

n_row <- nrow(metrics_tab)
row_pos <- tibble(
  id = metrics_tab$Method,
  row_i = seq_len(n_row),
  colour_background = seq_len(n_row) %% 2 == 1,
  ysep = row_space,
  y = -(seq_len(n_row) * row_height + cumsum(rep(row_space, n_row))),
)
row_pos$ymin <- row_pos$y - row_height / 2
row_pos$ymax <- row_pos$y + row_height / 2

column_info <- tibble(
  id = names(metrics_tab),
  group = c(
    "Text", "Image", "Text", "Text", "Score overall",
    rep("Removal of batch effects", 1 + n_metrics_batch),
    rep("Cell type label variance", 1 + n_metrics_bio)
  ),
  geom = c(
    "text", "image", "text", "text", "bar", "bar",
    rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)
  ),
  width = c(8, 2.5, 2, 1.5, 2, 2, rep(1, n_metrics_batch), 2, rep(1, n_metrics_bio)),
  overlay = FALSE
)

column_info$do_spacing <- c(FALSE, diff(as.integer(factor(column_info$group))) != 0)
column_info$xsep <- ifelse(column_info$do_spacing, col_bigspace, col_space)
column_info$xsep[1] <- 0
column_info$xmax <- cumsum(column_info$width + column_info$xsep)
column_info$xmin <- column_info$xmax - column_info$width
column_info$x <- column_info$xmin + column_info$width / 2

palettes <- list(
  `Score overall` = "YlGnBu",
  `Removal of batch effects` = "BuPu",
  `Cell type label variance` = "RdPu"
)

rank_fill <- function(values, palette_name) {
  n_ok <- sum(!is.na(values))
  pal <- colorRampPalette(rev(brewer.pal(9, palette_name)))(n_ok)
  pal[rank(values, ties.method = "average", na.last = "keep")]
}

# ---- 第七步：圆点、色条、文字、图标的绘图表 --------------------------------
ind_circle <- which(column_info$geom == "circle")
dat_circle <- as.matrix(metrics_tab[, ind_circle, drop = FALSE])
circle_data <- tibble(
  label = rep(colnames(dat_circle), each = n_row),
  x0 = rep(column_info$x[ind_circle], each = n_row),
  y0 = rep(row_pos$y, times = ncol(dat_circle)),
  value = as.vector(dat_circle),
  r = row_height / 2 * sqrt(as.vector(dat_circle))
)
circle_data <- circle_data %>%
  group_by(label) %>%
  mutate(r = rescale(r, to = c(0.05, 0.55), from = range(r, na.rm = TRUE))) %>%
  ungroup()
circle_colors <- unlist(lapply(seq_len(ncol(dat_circle)), function(i) {
  grp <- column_info$group[ind_circle[i]]
  rank_fill(dat_circle[, i], palettes[[grp]])
}))
circle_data$colors <- circle_colors
circle_data <- circle_data[!is.na(circle_data$value), ]

ind_bar <- which(column_info$geom == "bar")
dat_bar <- as.matrix(metrics_tab[, ind_bar, drop = FALSE])
rect_data <- tibble(
  label = rep(colnames(dat_bar), each = n_row),
  xmin = rep(column_info$xmin[ind_bar], each = n_row),
  xmax = rep(column_info$xmax[ind_bar], each = n_row),
  ymin = rep(row_pos$ymin, times = ncol(dat_bar)),
  ymax = rep(row_pos$ymax, times = ncol(dat_bar)),
  xwidth = rep(column_info$width[ind_bar], each = n_row),
  value = as.vector(dat_bar)
) %>%
  mutate(
    xmin = xmin + (1 - value) * xwidth * 0,
    xmax = xmax - (1 - value) * xwidth * 1
  )
rect_data$colors <- unlist(lapply(seq_len(ncol(dat_bar)), function(i) {
  grp <- column_info$group[ind_bar[i]]
  rank_fill(dat_bar[, i], palettes[[grp]])
}))
rect_data <- rect_data[!is.na(rect_data$value), ]

ind_text <- which(column_info$geom == "text")
dat_text <- as.matrix(metrics_tab[, ind_text, drop = FALSE])
text_data <- tibble(
  label_value = as.vector(dat_text),
  group = rep(colnames(dat_text), each = n_row),
  xmin = rep(column_info$xmin[ind_text], each = n_row),
  xmax = rep(column_info$xmax[ind_text], each = n_row),
  ymin = rep(row_pos$ymin, times = ncol(dat_text)),
  ymax = rep(row_pos$ymax, times = ncol(dat_text)),
  size = 4,
  fontface = "plain",
  colors = "black",
  angle = 0,
  hjust = 0.5,
  vjust = 0.5
)
text_data$colors[text_data$label_value == "HVG"] <- "darkgreen"
text_data$colors[text_data$label_value == "FULL"] <- "grey30"
text_data$label_value[text_data$label_value == "scaled"] <- "+"
text_data$label_value[text_data$label_value == "unscaled"] <- "-"
is_sign <- text_data$label_value %in% c("+", "-")
text_data$size[is_sign] <- 5
text_data$fontface[is_sign] <- "bold"

header_df <- column_info %>% filter(.data$id != "Method")
segment_data <- tibble(
  x = header_df$x,
  xend = header_df$x,
  y = -0.3,
  yend = -0.1,
  size = 0.5,
  colour = "black",
  linetype = "solid"
)
text_data <- bind_rows(
  text_data,
  tibble(
    xmin = header_df$x,
    xmax = header_df$x,
    ymin = 0,
    ymax = -0.5,
    angle = 30,
    vjust = 0,
    hjust = 0,
    label_value = header_df$id,
    size = 3,
    fontface = "plain",
    colors = "black",
    group = "header"
  )
)

# ---- 第八步：底部图例（Output / Scaling / Ranking / Score / 星号） ----------
minimum_x <- min(column_info$xmin, text_data$xmin, na.rm = TRUE)
maximum_x <- max(column_info$xmax, text_data$xmax, na.rm = TRUE)
minimum_y <- min(row_pos$ymin, text_data$ymin, na.rm = TRUE)
maximum_y <- max(row_pos$ymax, text_data$ymax, na.rm = TRUE)

leg_max_y <- minimum_y - 0.5
x_min_output <- minimum_x + 0.5
x_min_scaling <- minimum_x + 5.5
x_min_ranking <- minimum_x + 10.5
x_min_score <- minimum_x + 17

text_data <- bind_rows(
  text_data,
  tibble(
    xmin = x_min_output, xmax = x_min_output + 2,
    ymin = leg_max_y - 1, ymax = leg_max_y,
    label_value = "Output", hjust = 0, vjust = 0, fontface = "bold", size = 3,
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = x_min_output + 1.5, xmax = x_min_output + 3,
    ymin = c(leg_max_y - 2.2, leg_max_y - 3.4, leg_max_y - 4.6),
    ymax = c(leg_max_y - 2.2, leg_max_y - 3.4, leg_max_y - 4.6),
    label_value = c("gene", "embed", "graph"),
    hjust = 0, vjust = 0, fontface = "plain", size = 3,
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = x_min_scaling, xmax = x_min_scaling + 2,
    ymin = leg_max_y - 1, ymax = leg_max_y,
    label_value = "Scaling", hjust = 0, vjust = 0, fontface = "bold", size = 3,
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = c(x_min_scaling, x_min_scaling + 1, x_min_scaling, x_min_scaling + 1),
    xmax = c(x_min_scaling + 0.5, x_min_scaling + 3, x_min_scaling + 0.5, x_min_scaling + 3),
    ymin = c(leg_max_y - 2, leg_max_y - 2, leg_max_y - 3, leg_max_y - 3),
    ymax = c(leg_max_y - 1, leg_max_y - 1, leg_max_y - 2, leg_max_y - 2),
    label_value = c("+", ": scaled", "-", ": unscaled"),
    hjust = 0, vjust = 0,
    fontface = c("bold", "plain", "bold", "plain"),
    size = c(5, 3, 5, 3),
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = x_min_ranking, xmax = x_min_ranking + 2,
    ymin = leg_max_y - 1, ymax = leg_max_y,
    label_value = "Ranking", hjust = 0, vjust = 0, fontface = "bold", size = 3,
    colors = "black", angle = 0, group = "legend"
  )
)

rank_min_x <- c(
  `Score overall` = x_min_ranking,
  `Removal of batch effects` = x_min_ranking + 1,
  `Cell type label variance` = x_min_ranking + 2
)
rank_groups <- unique(column_info$group[column_info$geom == "bar"])
for (rg in rank_groups) {
  rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)
  rect_data <- bind_rows(
    rect_data,
    tibble(
      xmin = rank_min_x[[rg]],
      xmax = rank_min_x[[rg]] + 0.8,
      ymin = seq(leg_max_y - 4, leg_max_y - 2, by = 0.5),
      ymax = seq(leg_max_y - 3.5, leg_max_y - 1.5, by = 0.5),
      colors = rank_palette,
      label = "rank_legend",
      value = NA_real_,
      xwidth = 0.8
    )
  )
}
leg_max_x <- x_min_ranking + 2
arrow_data <- tibble(
  x = leg_max_x + 1.5, xend = leg_max_x + 1.5,
  y = leg_max_y - 4, yend = leg_max_y - 1.5,
  size = 0.5, colour = "black", linetype = "solid"
)
text_data <- bind_rows(
  text_data,
  tibble(
    xmin = leg_max_x + 2, xmax = leg_max_x + 2.5,
    ymin = c(leg_max_y - 2, leg_max_y - 4),
    ymax = c(leg_max_y - 1.5, leg_max_y - 3.5),
    label_value = c("1", as.character(n_row)),
    hjust = 0, vjust = 0, size = 2.5, fontface = "plain",
    colors = "black", angle = 0, group = "legend"
  )
)

cir_legend_dat <- tibble(value = seq(0, 1, by = 0.2), r = row_height / 2 * seq(0, 1, by = 0.2))
cir_legend_dat$r <- rescale(
  cir_legend_dat$r,
  to = c(0.05, 0.55),
  from = range(cir_legend_dat$r, na.rm = TRUE)
)
x0 <- numeric(nrow(cir_legend_dat))
cir_legend_space <- 0.1
for (i in seq_len(nrow(cir_legend_dat))) {
  if (i == 1) {
    x0[i] <- x_min_score + cir_legend_space + cir_legend_dat$r[i]
  } else {
    x0[i] <- x0[i - 1] + cir_legend_dat$r[i - 1] + cir_legend_space + cir_legend_dat$r[i]
  }
}
cir_legend_dat$x0 <- x0
cir_legend_min_y <- leg_max_y - 4
cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r
cir_legend_dat$colors <- NA_character_
cir_legend_dat$label <- "score_legend"

circle_data <- bind_rows(
  circle_data,
  cir_legend_dat %>% transmute(label, x0, y0, r, colors, value)
)
circle_data$colors[is.na(circle_data$colors)] <- "white"
text_data <- bind_rows(
  text_data,
  tibble(
    xmin = x_min_score, xmax = max(cir_legend_dat$x0),
    ymin = leg_max_y - 1, ymax = leg_max_y,
    label_value = "Score", hjust = 0, vjust = 0, fontface = "bold", size = 3,
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = cir_legend_dat$x0 - cir_legend_dat$r,
    xmax = cir_legend_dat$x0 + cir_legend_dat$r,
    ymin = cir_legend_min_y,
    ymax = cir_legend_min_y + 3,
    hjust = 0.5, vjust = 0, size = 2.5, fontface = "plain",
    label_value = ifelse(cir_legend_dat$value %in% c(0, 1), paste0(cir_legend_dat$value * 100, "%"), ""),
    colors = "black", angle = 0, group = "legend"
  ),
  tibble(
    xmin = x_min_output,
    xmax = x_min_output + 8,
    ymin = leg_max_y - 7.2,
    ymax = leg_max_y - 5.4,
    label_value = "* uses coarse cell type labels as input",
    hjust = 0, vjust = 1, fontface = "plain", size = 3,
    colors = "black", angle = 0, group = "legend"
  )
)

# ---- 第九步：Output 图标（gene 矩阵 / embed 散点 / graph 网络） -------------
# 官方脚本用 cowplot::draw_image 贴 PNG。当前环境没有 png/ggimage，
# 这里按官方图标的几何用 ggplot 图层重画，放在同一 1×1 数据坐标盒里。
icon_x <- column_info$x[column_info$geom == "image"]
icon_rows <- tibble(
  output = metrics_tab$Output,
  x = icon_x,
  y = row_pos$y
)
legend_icons <- tibble(
  output = c("gene", "embed", "graph"),
  x = x_min_output + 0.5,
  y = c(leg_max_y - 2, leg_max_y - 3.2, leg_max_y - 4.4)
)
all_icons <- bind_rows(icon_rows, legend_icons)

gene_icons <- all_icons %>% filter(output == "gene")
embed_icons <- all_icons %>% filter(output == "embed")
graph_icons <- all_icons %>% filter(output == "graph")

# 4×4 白格 + 黑间隔，比例接近官方 matrix.png
# 官方 matrix.png 是黑间隔 + 透明方格；白底上看起来就是 4×4 白格。
gene_outer <- 0.42
gene_cell <- 0.155
gene_gap <- 0.053
gene_origin <- -2 * gene_cell - 1.5 * gene_gap
gene_tiles <- expand.grid(i = 0:3, j = 0:3)
gene_bg <- tibble(
  xmin = gene_icons$x - gene_outer,
  xmax = gene_icons$x + gene_outer,
  ymin = gene_icons$y - gene_outer,
  ymax = gene_icons$y + gene_outer
)
gene_rect <- bind_rows(lapply(seq_len(nrow(gene_icons)), function(k) {
  tibble(
    xmin = gene_icons$x[k] + gene_origin + gene_tiles$i * (gene_cell + gene_gap),
    xmax = gene_icons$x[k] + gene_origin + gene_tiles$i * (gene_cell + gene_gap) + gene_cell,
    ymin = gene_icons$y[k] + gene_origin + gene_tiles$j * (gene_cell + gene_gap),
    ymax = gene_icons$y[k] + gene_origin + gene_tiles$j * (gene_cell + gene_gap) + gene_cell
  )
}))

embed_pts_unit <- tibble(
  px = c(0.18, 0.42, 0.68, 0.22, 0.50, 0.78, 0.28, 0.58, 0.80, 0.20, 0.40, 0.62, 0.24, 0.48, 0.76),
  py = c(0.22, 0.18, 0.26, 0.38, 0.34, 0.42, 0.52, 0.48, 0.58, 0.66, 0.62, 0.70, 0.82, 0.78, 0.86)
)
embed_points <- bind_rows(lapply(seq_len(nrow(embed_icons)), function(k) {
  tibble(
    x = embed_icons$x[k] - 0.32 + 0.68 * embed_pts_unit$px,
    y = embed_icons$y[k] - 0.32 + 0.68 * embed_pts_unit$py
  )
}))
embed_axes <- bind_rows(lapply(seq_len(nrow(embed_icons)), function(k) {
  tibble(
    x = c(embed_icons$x[k] - 0.32, embed_icons$x[k] - 0.32),
    xend = c(embed_icons$x[k] - 0.32, embed_icons$x[k] + 0.36),
    y = c(embed_icons$y[k] - 0.32, embed_icons$y[k] - 0.32),
    yend = c(embed_icons$y[k] + 0.36, embed_icons$y[k] - 0.32)
  )
}))

graph_nodes_unit <- tibble(
  nx = c(0.126, 0.690, 0.459, 0.872, 0.184, 0.874),
  ny = 1 - c(0.192, 0.125, 0.461, 0.438, 0.829, 0.875)
)
graph_edges_idx <- tibble(a = c(1, 1, 1, 2, 3, 3, 4, 5), b = c(2, 3, 5, 3, 4, 6, 5, 6))
graph_edges <- bind_rows(lapply(seq_len(nrow(graph_icons)), function(k) {
  tibble(
    x = graph_icons$x[k] - 0.40 + 0.80 * graph_nodes_unit$nx[graph_edges_idx$a],
    y = graph_icons$y[k] - 0.40 + 0.80 * graph_nodes_unit$ny[graph_edges_idx$a],
    xend = graph_icons$x[k] - 0.40 + 0.80 * graph_nodes_unit$nx[graph_edges_idx$b],
    yend = graph_icons$y[k] - 0.40 + 0.80 * graph_nodes_unit$ny[graph_edges_idx$b]
  )
}))
graph_nodes <- bind_rows(lapply(seq_len(nrow(graph_icons)), function(k) {
  tibble(
    x0 = graph_icons$x[k] - 0.40 + 0.80 * graph_nodes_unit$nx,
    y0 = graph_icons$y[k] - 0.40 + 0.80 * graph_nodes_unit$ny,
    r = 0.075
  )
}))

# ---- 第十步：拼 ggplot ------------------------------------------------------
text_data <- text_data %>%
  mutate(
    angle2 = angle / 360 * 2 * pi,
    cosa = round(cos(angle2), 2),
    sina = round(sin(angle2), 2),
    alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) +
      ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
    alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) +
      ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
    x = (1 - alphax) * xmin + alphax * xmax,
    y = (1 - alphay) * ymin + alphay * ymax
  ) %>%
  filter(label_value != "")

text_left <- text_data %>% filter(group == "Method")
text_left$x <- text_left$x - 3
text_other <- text_data %>% filter(group != "Method")

bg_rows <- row_pos %>% filter(colour_background)

g <- ggplot() +
  coord_equal(expand = FALSE) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_size_identity() +
  scale_linewidth_identity() +
  scale_linetype_identity() +
  theme_void(base_family = plot_font_family) +
  theme(plot.margin = margin(0, 0, 0, 0))

if (nrow(bg_rows) > 0) {
  g <- g + geom_rect(
    data = bg_rows,
    aes(
      xmin = min(column_info$xmin) - 0.25,
      xmax = max(column_info$xmax) + 0.25,
      ymin = ymin - row_space / 2,
      ymax = ymax + row_space / 2
    ),
    fill = "#DDDDDD",
    colour = NA
  )
}

g <- g +
  geom_circle(
    data = circle_data,
    aes(x0 = x0, y0 = y0, r = r, fill = colors),
    colour = "black",
    linewidth = 0.25,
    na.rm = TRUE
  ) +
  geom_rect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors),
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = text_other,
    aes(
      x = x, y = y, label = label_value, colour = colors,
      hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle
    ),
    family = plot_font_family
  ) +
  geom_text(
    data = text_left,
    aes(
      x = x, y = y, label = label_value, colour = colors,
      vjust = vjust, size = size, fontface = fontface, angle = angle
    ),
    hjust = 0,
    family = plot_font_family
  ) +
  geom_segment(
    data = segment_data,
    aes(x = x, xend = xend, y = y, yend = yend, linewidth = size, colour = colour, linetype = linetype)
  ) +
  geom_segment(
    data = arrow_data,
    aes(x = x, xend = xend, y = y, yend = yend, linewidth = size, colour = colour, linetype = linetype),
    arrow = arrow(length = unit(0.1, "cm")),
    lineend = "round",
    linejoin = "bevel"
  ) +
  geom_rect(
    data = gene_bg,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "black",
    colour = NA
  ) +
  geom_rect(
    data = gene_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    colour = NA
  ) +
  geom_segment(
    data = embed_axes,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = "black",
    linewidth = 0.28,
    arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
    lineend = "round"
  ) +
  geom_point(
    data = embed_points,
    aes(x = x, y = y),
    colour = "black",
    size = 0.35
  ) +
  geom_segment(
    data = graph_edges,
    aes(x = x, xend = xend, y = y, yend = yend),
    colour = "black",
    linewidth = 0.32,
    lineend = "round"
  ) +
  geom_circle(
    data = graph_nodes,
    aes(x0 = x0, y0 = y0, r = r),
    fill = "black",
    colour = "black",
    linewidth = 0
  )

minimum_x <- minimum_x - 2
maximum_x <- maximum_x + 5
minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE)) - 2
# 官方 knit_table：maximum_y + 4，给 30° 表头留空。
# 再加一段给 “SUPPLEMENTARY FIGURES”，标题跟方法名左对齐。
# A3 页的右侧留白是官方 ggsave(297×420 mm) 就有的，不是再额外撑开。
heading_x <- min(text_left$x, na.rm = TRUE)
maximum_y <- maximum_y + 4 + 1.6
g <- g +
  expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y)) +
  annotate(
    "text",
    x = heading_x,
    y = maximum_y - 0.2,
    label = "SUPPLEMENTARY FIGURES",
    hjust = 0,
    vjust = 1,
    size = 5.4,
    fontface = "bold",
    family = plot_font_family,
    colour = "black"
  ) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ---- 第十一步：导出 PNG / PDF ----------------------------------------------
# 官方 ggsave 用 A3（297×420 mm），右侧会留白，和期刊补图页一致。
png_a3 <- file.path(out_dir, "hlca_scib_supplementary_figure1.png")
pdf_a3 <- file.path(out_dir, "hlca_scib_supplementary_figure1.pdf")
png_compact <- file.path(out_dir, "hlca_scib_supplementary_figure1_compact.png")

ggsave(
  png_a3,
  plot = g,
  device = ragg::agg_png,
  width = 297,
  height = 420,
  units = "mm",
  dpi = 300,
  background = "white",
  limitsize = FALSE
)
ggsave(
  pdf_a3,
  plot = g,
  device = cairo_pdf,
  width = 297,
  height = 420,
  units = "mm",
  bg = "white",
  limitsize = FALSE
)

unit_mm <- 5.2
ggsave(
  png_compact,
  plot = g,
  device = ragg::agg_png,
  width = (maximum_x - minimum_x) * unit_mm,
  height = (maximum_y - minimum_y) * unit_mm,
  units = "mm",
  dpi = 300,
  background = "white",
  limitsize = FALSE
)

file.copy(png_compact, file.path(draft_dir, "preview.png"), overwrite = TRUE)

message("Wrote:")
message("  ", summary_out)
message("  ", png_a3)
message("  ", pdf_a3)
message("  ", png_compact)
message("n methods plotted: ", n_row)

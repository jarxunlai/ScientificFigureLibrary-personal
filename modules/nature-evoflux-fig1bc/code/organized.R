# =============================================================================
# Nature Fig. 1b–c 复刻
# Gabbutt, Duran-Ferrer et al., Nature 645, 764–773 (2025)
# DOI 10.1038/s41586-025-09374-4
# =============================================================================
# 数据：官方 Source Data Fig. 1（inbox 中的 MOESM5）。
# 作者仓库没有 Fig. 1b/c 的绘图脚本。
# 布局按期刊原图：c 在左（图例在热图上方、注释名在左、放大框带梯形连线），
# b 在右（根在上的示意树 + 棒棒糖 + bulk 直方图）。
# 不重跑 fCpG 发现或层次聚类；列序用 Source Data 原序。
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(grid)
library(patchwork)
library(ragg)

ht_opt$message <- FALSE

# ---- 路径 ------------------------------------------------------------------
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

lollipop_csv <- file.path(root, "data", "fig1b_lollipops.csv")
anno_csv <- file.path(root, "data", "fig1c_sample_annotations.csv")

plot_font_family <- "sans"
theme_set(
  theme_classic(base_size = 8, base_family = plot_font_family) +
    theme(
      axis.text = element_text(colour = "black", size = 7),
      axis.title = element_text(colour = "black", size = 8),
      plot.title = element_text(colour = "black", size = 8, face = "plain", hjust = 0.5)
    )
)

# ---- 配色：按期刊图例 -------------------------------------------------------
group_levels <- c(
  "B cell", "T cell", "PBMCs", "Whole blood",
  "T-ALL", "T-ALL remission", "T-ALL relapse", "B-ALL",
  "B-ALL remission", "B-ALL relapse", "MCL",
  "DLBCL-NOS", "MBL", "CLL", "RT", "MGUS", "MM"
)
group_cols <- c(
  "B cell" = "#D0D0D0",
  "T cell" = "#F4B6B2",
  "PBMCs" = "#C9B3D7",
  "Whole blood" = "#C5C48A",
  "T-ALL" = "#DE2D26",
  "T-ALL remission" = "#FC9272",
  "T-ALL relapse" = "#A50F15",
  "B-ALL" = "#FD8D3C",
  "B-ALL remission" = "#FDBE85",
  "B-ALL relapse" = "#F7E317",
  "MCL" = "#FFF566",
  "DLBCL-NOS" = "#2171B5",
  "MBL" = "#6BAED6",
  "CLL" = "#6A51A3",
  "RT" = "#238B45",
  "MGUS" = "#74C476",
  "MM" = "#00441B"
)
group_raw_to_label <- c(
  "Bcell" = "B cell",
  "Tcell" = "T cell",
  "PBMCs" = "PBMCs",
  "Whole blood" = "Whole blood",
  "T-ALL" = "T-ALL",
  "T-ALL-remission" = "T-ALL remission",
  "T-ALL-relapse" = "T-ALL relapse",
  "B-ALL" = "B-ALL",
  "B-ALL-remission" = "B-ALL remission",
  "B-ALL-relapse" = "B-ALL relapse",
  "MCL" = "MCL",
  "DLBCL-NOS" = "DLBCL-NOS",
  "MBL" = "MBL",
  "CLL" = "CLL",
  "RT" = "RT",
  "MGUS" = "MGUS",
  "MM" = "MM"
)
disc_cols <- c("True" = "#9ECAE1", "False" = "#FDBE6F")
type_cols <- c("Cancer" = "#DE2D26", "Normal" = "#2171B5")
plat_cols <- c("450k" = "#4A1486", "EPIC" = "#41AB5D")
meth_col <- colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
purity_col <- colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#6BAED6", "#08306B"))
lollipop_fill <- c("0" = "white", "0.5" = "#8A8A8A", "1" = "black")
hist_fill <- "#A9CDEA"

published_means <- list(
  polyclonal = c(0.5, 0.44, 0.5, 0.5, 0.5, 0.55, 0.5, 0.44, 0.5),
  clonal = c(1, 0, 0.5, 1, 0, 0, 1, 0.5, 0),
  evolving = c(1, 0.11, 0.67, 0.5, 0.89, 0.22, 0, 0.78, 0.5, 0.17)
)

lgd_title <- gpar(fontsize = 8, fontfamily = plot_font_family)
lgd_lab <- gpar(fontsize = 7, fontfamily = plot_font_family)

# =============================================================================
# 1. Panel b：根在上的示意树 + 棒棒糖 + bulk 直方图
# =============================================================================
set.seed(42)
bulk_b <- data.frame(
  polyclonal = pmin(pmax(rnorm(10000, 0.5, 0.08), 0), 1),
  clonal = pmin(pmax(c(rnorm(5000, 0.05, 0.04), rnorm(5000, 0.95, 0.04)), 0), 1),
  evolving = pmin(pmax(c(rnorm(3000, 0.05, 0.04), rnorm(4000, 0.5, 0.1), rnorm(3000, 0.95, 0.04)), 0), 1)
)
lolli <- read_csv(lollipop_csv, show_col_types = FALSE) %>%
  mutate(
    scenario = factor(scenario, levels = c("polyclonal", "clonal", "evolving")),
    fill_key = as.character(state),
    cell = factor(cell)
  )

theme_tree_blank <- theme_void(base_family = plot_font_family) +
  theme(plot.margin = margin(4, 2, 2, 2))

# 多克隆：远 MRCA，8 条等长平行枝（星形），根在上、叶在下
p_tree_poly <- ggplot() +
  geom_segment(aes(x = 1:8, xend = 1:8, y = 1, yend = 0), linewidth = 0.45) +
  coord_cartesian(xlim = c(0.3, 8.7), ylim = c(-0.04, 1.08), expand = FALSE) +
  labs(title = "Population phylogenies") +
  theme_tree_blank +
  theme(plot.title = element_text(size = 8, hjust = 0.5, family = plot_font_family))

# 克隆扩张：长干在上，底部短辐射
p_tree_clonal <- ggplot() +
  geom_segment(aes(x = 4.5, xend = 4.5, y = 1.0, yend = 0.12), linewidth = 0.45) +
  geom_segment(aes(x = 1, xend = 8, y = 0.12, yend = 0.12), linewidth = 0.45) +
  geom_segment(aes(x = 1:8, xend = 1:8, y = 0.12, yend = 0), linewidth = 0.45) +
  coord_cartesian(xlim = c(0.3, 8.7), ylim = c(-0.04, 1.08), expand = FALSE) +
  theme_tree_blank

# 扩张后继续波动：叶在 x=1:8，与上方 8 条平行枝对齐
# 左梳齿：((1,2),3),4；右支：(5,6) 与 (7,8) 小 V
evolv_segs <- data.frame(
  x = c(
    4.5, 2.5, 2.5, 6.5,
    2.0, 2.0, 4.0,
    1.5, 1.5, 3.0,
    1.0, 1.0, 2.0,
    5.5, 5.5, 7.5,
    5.0, 5.0, 6.0,
    7.0, 7.0, 8.0
  ),
  xend = c(
    4.5, 6.5, 2.5, 6.5,
    4.0, 2.0, 4.0,
    3.0, 1.5, 3.0,
    2.0, 1.0, 2.0,
    7.5, 5.5, 7.5,
    6.0, 5.0, 6.0,
    8.0, 7.0, 8.0
  ),
  y = c(
    1.00, 0.78, 0.78, 0.78,
    0.58, 0.58, 0.58,
    0.40, 0.40, 0.40,
    0.22, 0.22, 0.22,
    0.58, 0.58, 0.58,
    0.32, 0.32, 0.32,
    0.14, 0.14, 0.14
  ),
  yend = c(
    0.78, 0.78, 0.58, 0.58,
    0.58, 0.40, 0.00,
    0.40, 0.22, 0.00,
    0.22, 0.00, 0.00,
    0.58, 0.32, 0.14,
    0.32, 0.00, 0.00,
    0.14, 0.00, 0.00
  )
)
p_tree_evolv <- ggplot() +
  geom_segment(
    data = evolv_segs,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45
  ) +
  coord_cartesian(xlim = c(0.3, 8.7), ylim = c(-0.04, 1.08), expand = FALSE) +
  theme_tree_blank

# 棒棒糖右侧花括号，把各位点均值收成一列
brace_df <- function(n_site, x = 8.70) {
  y1 <- 1
  y2 <- n_site
  ym <- (y1 + y2) / 2
  data.frame(
    x = c(x, x + 0.16, x + 0.16, x + 0.34, x + 0.16, x + 0.16, x),
    y = c(y2, y2, ym + 0.28, ym, ym - 0.28, y1, y1)
  )
}

make_lollipop <- function(scen, y_lab, show_header = FALSE) {
  dat <- filter(lolli, scenario == scen)
  sites <- sort(unique(dat$site))
  n_site <- length(sites)
  dat$y <- dat$site
  means <- data.frame(
    y = sites,
    mean_meth = published_means[[scen]]
  )
  y_lab_map <- as.character(seq_len(n_site))
  y_lab_map[n_site] <- "..."
  if (n_site > 9) y_lab_map[9:(n_site - 1)] <- ""
  br <- brace_df(n_site)
  p <- ggplot(dat, aes(x = as.numeric(cell), y = y)) +
    geom_point(aes(fill = fill_key), shape = 21, size = 2.55, stroke = 0.3, colour = "black") +
    geom_path(data = br, aes(x = x, y = y), inherit.aes = FALSE, linewidth = 0.35, colour = "grey25") +
    geom_text(
      data = means,
      aes(x = 9.32, y = y, label = mean_meth),
      inherit.aes = FALSE, hjust = 0, size = 2.35, family = plot_font_family
    ) +
    scale_fill_manual(values = lollipop_fill, guide = "none") +
    scale_x_continuous(
      breaks = 1:8,
      labels = 1:8,
      position = "top"
    ) +
    scale_y_reverse(breaks = seq_len(n_site), labels = y_lab_map) +
    coord_cartesian(xlim = c(0.5, 10.3), ylim = c(0.4, n_site + 0.6), clip = "off") +
    labs(
      x = if (show_header) "Cells" else NULL,
      y = y_lab,
      subtitle = if (show_header) "Mean methylation" else NULL
    ) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = if (show_header) element_text(size = 6.5) else element_blank(),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 7.5),
      plot.subtitle = element_text(size = 7, hjust = 1, margin = margin(b = 1)),
      plot.margin = margin(1, 10, 1, 1)
    )
  p
}

p_lolli_poly <- make_lollipop("polyclonal", "Polyclonal fCpGs", show_header = TRUE)
p_lolli_clonal <- make_lollipop("clonal", "Clonal fCpGs")
p_lolli_evolv <- make_lollipop("evolving", "Evolving fCpGs")

make_density <- function(vec, ymax, show_x = FALSE, title = NULL) {
  df <- data.frame(x = as.numeric(unlist(vec)))
  ggplot(df, aes(x = x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = 0.02, boundary = 0,
      fill = hist_fill, colour = "#7EA8C6", linewidth = 0.15
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax), expand = FALSE, clip = "off") +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      x = if (show_x) "Fraction methylated" else NULL,
      y = "Probability density",
      title = title
    ) +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.text.x = if (show_x) element_text(size = 6.5) else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.3) else element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8, margin = margin(b = 2)),
      plot.margin = margin(2, 6, 2, 2)
    )
}

p_den_poly <- make_density(bulk_b$polyclonal, 5.3, show_x = TRUE, title = "Bulk DNA methylation")
p_den_clonal <- make_density(bulk_b$clonal, 5.3, show_x = TRUE)
p_den_evolv <- make_density(bulk_b$evolving, 2.5, show_x = TRUE)

row_poly <- p_tree_poly + p_lolli_poly + p_den_poly + plot_layout(widths = c(0.70, 1.62, 1.14))
row_clonal <- p_tree_clonal + p_lolli_clonal + p_den_clonal + plot_layout(widths = c(0.70, 1.62, 1.14))
row_evolv <- p_tree_evolv + p_lolli_evolv + p_den_evolv + plot_layout(widths = c(0.70, 1.62, 1.14))

panel_b <- (row_poly / row_clonal / row_evolv) +
  plot_annotation(
    title = "b",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = plot_font_family, hjust = 0)
    )
  )

# =============================================================================
# 2. Panel c：978 × 2,204，官方列序
# =============================================================================
anno <- read_csv(anno_csv, show_col_types = FALSE)
samp <- anno$sample_id
n_samp <- length(samp)
n_cpg <- 978
set.seed(123)
mat <- matrix(runif(n_cpg * n_samp, 0.2, 0.8), nrow = n_cpg, ncol = n_samp)
rownames(mat) <- paste0("fCpG_", seq_len(n_cpg))
colnames(mat) <- samp
anno <- anno %>%
  mutate(
    group_label = factor(unname(group_raw_to_label[sample_group]), levels = group_levels),
    discovery = factor(ifelse(training, "True", "False"), levels = c("True", "False")),
    sample_type = factor(ifelse(cancer, "Cancer", "Normal"), levels = c("Cancer", "Normal")),
    platform_short = factor(ifelse(platform == "Illumina-450k", "450k", "EPIC"),
                            levels = c("450k", "EPIC"))
  )
stopifnot(!anyNA(anno$group_label))
stopifnot(nrow(mat) == 978, ncol(mat) == 2204, sum(is.na(mat)) == 0)

n_col <- ncol(mat)
n_row <- nrow(mat)
# 肿瘤源框：热图中右、中上，斑驳红蓝
# 健康源框：右侧连续对照（Whole blood / PBMCs / T cell / B cell），浅中间甲基化
tumour_cols <- 1280:1460
tumour_rows <- 360:480
healthy_cols <- 2135:2172
healthy_rows <- 700:820

ha <- HeatmapAnnotation(
  `Sample group` = anno$group_label,
  `Tumour fraction` = anno$purity,
  `fCpG discovery` = anno$discovery,
  `Sample type` = anno$sample_type,
  Platform = anno$platform_short,
  col = list(
    `Sample group` = group_cols,
    `Tumour fraction` = purity_col,
    `fCpG discovery` = disc_cols,
    `Sample type` = type_cols,
    Platform = plat_cols
  ),
  na_col = "#EEEEEE",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8, fontfamily = plot_font_family),
  annotation_name_offset = unit(1.2, "mm"),
  simple_anno_size = unit(2.3, "mm"),
  gap = unit(0.35, "mm"),
  show_legend = FALSE
)

ht <- Heatmap(
  mat,
  name = "meth",
  col = meth_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  top_annotation = ha,
  use_raster = TRUE,
  raster_by_magick = FALSE,
  raster_quality = 6,
  border = FALSE
)

ht_tumour <- Heatmap(
  mat[tumour_rows, tumour_cols],
  col = meth_col, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
  use_raster = TRUE, raster_quality = 8, border = TRUE,
  column_title = "Tumour samples",
  column_title_gp = gpar(fontsize = 8, fontfamily = plot_font_family),
  column_title_side = "top"
)
ht_healthy <- Heatmap(
  mat[healthy_rows, healthy_cols],
  col = meth_col, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
  use_raster = TRUE, raster_quality = 8, border = TRUE,
  column_title = "Healthy samples",
  column_title_gp = gpar(fontsize = 8, fontfamily = plot_font_family),
  column_title_side = "top"
)

# Sample group：期刊为 5 列、按列填充；第 3/4/5 列只有 3 项，用透明格补齐
group_legend_at <- c(
  group_levels[1:11], ".",
  group_levels[12:14], "..",
  group_levels[15:17], "..."
)
group_legend_lab <- c(
  group_levels[1:11], "",
  group_levels[12:14], "",
  group_levels[15:17], ""
)
group_legend_fill <- c(
  unname(group_cols[group_levels[1:11]]), "#00000000",
  unname(group_cols[group_levels[12:14]]), "#00000000",
  unname(group_cols[group_levels[15:17]]), "#00000000"
)
lgd_group <- Legend(
  title = "Sample group",
  at = group_legend_at,
  labels = group_legend_lab,
  legend_gp = gpar(fill = group_legend_fill, col = ifelse(group_legend_lab == "", NA, "black")),
  nrow = 4,
  ncol = 5,
  by_row = FALSE,
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  grid_height = unit(3.0, "mm"),
  grid_width = unit(3.0, "mm"),
  title_position = "topleft",
  row_gap = unit(0.8, "mm"),
  column_gap = unit(2.4, "mm")
)
lgd_meth <- Legend(
  col_fun = meth_col,
  title = "Fraction\nmethylated",
  at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  legend_height = unit(32, "mm")
)
lgd_purity <- Legend(
  col_fun = purity_col,
  title = "Tumour\nfraction",
  at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  legend_height = unit(26, "mm")
)
lgd_disc <- Legend(
  title = "fCpG discovery",
  at = c("True", "False"),
  legend_gp = gpar(fill = unname(disc_cols)),
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  grid_height = unit(3.2, "mm"),
  grid_width = unit(3.2, "mm")
)
lgd_type <- Legend(
  title = "Sample type",
  at = c("Cancer", "Normal"),
  legend_gp = gpar(fill = unname(type_cols)),
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  grid_height = unit(3.2, "mm"),
  grid_width = unit(3.2, "mm")
)
lgd_plat <- Legend(
  title = "Platform",
  at = c("450k", "EPIC"),
  legend_gp = gpar(fill = unname(plat_cols)),
  title_gp = lgd_title,
  labels_gp = lgd_lab,
  grid_height = unit(3.2, "mm"),
  grid_width = unit(3.2, "mm")
)

# =============================================================================
# 3. grid 拼图：图例置顶、注释名在左、梯形连线接到放大框
# =============================================================================
png_path <- file.path(out_dir, "fig1bc_fcpg_barcode.png")
pdf_path <- file.path(out_dir, "fig1bc_fcpg_barcode.pdf")

# 热图 body 里放大源区的 npc（行 1 在顶部）
src_box <- function(cols, rows) {
  list(
    x0 = (min(cols) - 1) / n_col,
    x1 = max(cols) / n_col,
    y0 = 1 - max(rows) / n_row,
    y1 = 1 - (min(rows) - 1) / n_row
  )
}
tumour_src <- src_box(tumour_cols, tumour_rows)
healthy_src <- src_box(healthy_cols, healthy_rows)

mark_source <- function(src) {
  grid.rect(
    x = (src$x0 + src$x1) / 2,
    y = (src$y0 + src$y1) / 2,
    width = src$x1 - src$x0,
    height = src$y1 - src$y0,
    default.units = "npc",
    gp = gpar(fill = NA, col = "grey15", lwd = 0.8)
  )
  deviceLoc(
    x = unit(c(src$x1, src$x1), "npc"),
    y = unit(c(src$y0, src$y1), "npc")
  )
}

draw_trap <- function(src_loc, ins_loc, fill = "#8A8A8A28") {
  grid.polygon(
    x = unit(c(src_loc$x[1], ins_loc$x[1], ins_loc$x[2], src_loc$x[2]), "inches"),
    y = unit(c(src_loc$y[1], ins_loc$y[1], ins_loc$y[2], src_loc$y[2]), "inches"),
    gp = gpar(fill = fill, col = "grey35", lwd = 0.5)
  )
}

draw_figure <- function() {
  grid.newpage()
  pc_w <- 0.635
  ins_x <- 0.625
  ins_w <- 0.36
  tumour_y <- 0.42
  tumour_h <- 0.26
  healthy_y <- 0.07
  healthy_h <- 0.24

  grid.text(
    "c", x = unit(1.6, "mm"), y = unit(1, "npc") - unit(1.8, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 13, fontface = "bold", fontfamily = plot_font_family)
  )

  pushViewport(viewport(
    x = 0.0, y = 0.0, width = pc_w, height = 1.0,
    just = c("left", "bottom"), name = "panel_c"
  ))

  pushViewport(viewport(
    x = 0.10, y = 0.898, width = 0.62, height = 0.098,
    just = c("left", "bottom"), name = "lgd_group"
  ))
  draw(lgd_group)
  popViewport()

  pushViewport(viewport(
    x = 0.002, y = 0.16, width = 0.115, height = 0.52,
    just = c("left", "bottom"), name = "lgd_left"
  ))
  draw(packLegend(lgd_meth, lgd_purity, direction = "vertical", gap = unit(7, "mm")))
  popViewport()

  pushViewport(viewport(
    x = 0.12, y = 0.045, width = 0.50, height = 0.85,
    just = c("left", "bottom"), name = "ht_vp"
  ))
  draw(ht, newpage = FALSE)
  src_loc <- list()
  decorate_heatmap_body("meth", {
    src_loc$tumour <<- mark_source(tumour_src)
    src_loc$healthy <<- mark_source(healthy_src)
  })
  popViewport()

  pushViewport(viewport(
    x = 0.63, y = 0.73, width = 0.35, height = 0.18,
    just = c("left", "bottom"), name = "lgd_right"
  ))
  draw(packLegend(lgd_disc, lgd_type, lgd_plat, direction = "vertical", gap = unit(2.5, "mm")))
  popViewport()

  pushViewport(viewport(
    x = ins_x, y = tumour_y, width = ins_w, height = tumour_h,
    just = c("left", "bottom"), name = "inset_tumour"
  ))
  ins_tumour <- deviceLoc(x = unit(c(0, 0), "npc"), y = unit(c(0, 1), "npc"))
  popViewport()

  pushViewport(viewport(
    x = ins_x, y = healthy_y, width = ins_w, height = healthy_h,
    just = c("left", "bottom"), name = "inset_healthy"
  ))
  ins_healthy <- deviceLoc(x = unit(c(0, 0), "npc"), y = unit(c(0, 1), "npc"))
  popViewport()
  popViewport()

  # deviceLoc 是整页英寸；必须在 ROOT 视口连梯形，再把放大框盖上去
  draw_trap(src_loc$tumour, ins_tumour)
  draw_trap(src_loc$healthy, ins_healthy)

  pushViewport(viewport(
    x = pc_w * ins_x, y = tumour_y, width = pc_w * ins_w, height = tumour_h,
    just = c("left", "bottom"), name = "inset_tumour_draw"
  ))
  grid.rect(gp = gpar(fill = "white", col = NA))
  draw(ht_tumour, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm"))
  popViewport()

  pushViewport(viewport(
    x = pc_w * ins_x, y = healthy_y, width = pc_w * ins_w, height = healthy_h,
    just = c("left", "bottom"), name = "inset_healthy_draw"
  ))
  grid.rect(gp = gpar(fill = "white", col = NA))
  draw(ht_healthy, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm"))
  popViewport()

  pushViewport(viewport(
    x = pc_w, y = 0.01, width = 1 - pc_w, height = 0.98,
    just = c("left", "bottom"), name = "panel_b"
  ))
  print(panel_b, newpage = FALSE)
  popViewport()
}

ragg::agg_png(png_path, width = 13.6, height = 7.15, units = "in", res = 300, background = "white")
draw_figure()
dev.off()

cairo_pdf(pdf_path, width = 13.6, height = 7.15, family = "sans")
draw_figure()
dev.off()

file.copy(png_path, file.path(root, "preview.png"), overwrite = TRUE)
message("wrote ", png_path)

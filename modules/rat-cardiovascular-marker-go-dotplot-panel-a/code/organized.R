# Supplementary Figure 2A — 仅视觉复刻；输入数字为合成拟态，不复现上游分析。
# 从本模块目录执行: Rscript --vanilla code/plot.R

root <- "."
lab <- read.csv(file.path(root, "data/panel_a_labels.csv"), check.names = FALSE)
dots <- read.csv(file.path(root, "data/dotplot_synthetic.csv"))
nes <- read.csv(file.path(root, "data/go_nes_synthetic.csv"))

# 一、行列和来源断言：CSV 每行的 5 个 marker 必须连续向下排列。
cell <- lab$cell_type
genes <- as.vector(t(as.matrix(lab[paste0("gene_", 1:5)])))
owner <- rep(cell, each = 5)
stopifnot(nrow(lab) == 26L, length(genes) == 130L, !anyDuplicated(genes),
          identical(genes[1:5], c("Gsn", "Dcn", "Col3a1", "Opcml", "C7")),
          identical(tail(genes, 5), c("Cntn5", "Snhg11", "Stmn2", "Snap25", "Nefl")),
          nrow(dots) == 3380L, nrow(nes) == 78L,
          all(dots$gene %in% genes), all(dots$cell_type %in% cell),
          identical(nes$cell_type, rep(cell, each = 3)),
          identical(nes$rank, rep(1:3, 26)),
          all(dots$origin == "synthetic_visual_mimic"),
          all(nes$origin == "synthetic_visual_mimic"))
gene_index <- match(dots$gene, genes)
cell_index <- match(dots$cell_type, cell)
stopifnot(sum(dots$pct_expressing > 0 & dots$cell_type == owner[gene_index]) == 130L)

# 二、原图各轨道的固定宽度（单位 mm）：图例、marker 框、Cell type、NES 框、框外 GO 文字。
W <- 249; H <- 350
y_bottom <- 10; y_top <- 310; pitch <- (y_top - y_bottom) / length(genes)
dot_l <- 54; dot_r <- 107
nes_l <- 137; nes_r <- 155; go_x <- 158
x_dot <- dot_l + ((seq_along(cell) - 0.5) / length(cell)) * (dot_r - dot_l)
y_gene <- y_top - ((seq_along(genes) - 0.5) * pitch)
y_go <- y_top - (((rep(seq_along(cell), each = 3) - 1) * 5 +
                     rep(2:4, length(cell)) - 0.5) * pitch)
stopifnot(identical(owner[1:5], rep(cell[1], 5)),
          all(diff(y_gene) < 0), all(diff(y_go) < 0),
          y_gene[1] > y_go[1], y_gene[130] < y_go[78])

# 三、绘图。点大小/颜色与 NES 均为视觉拟态；标签与顺序取自原图。
draw_panel <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = "sans", ps = 9)
  plot.new()
  plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")
  cols <- grDevices::colorRampPalette(c("#E5EFF8", "#81AFD2", "#143F70"))(101)

  # 独立发布版本不保留原多面板的 a 标记。
  text((dot_l + dot_r)/2, 342, "Marker genes", cex = 1.12, font = 2)
  text(110, 342, "Cell type", cex = 1.1, font = 2, adj = 0)
  text(nes_l + 4, 342, "GO terms", cex = 1.1, font = 2, adj = 0)

  # 原图顶部旋转的横轴坐标文字为黑色；只有中间 Cell type 注释着色。
  for (i in seq_along(genes))
    text(dot_l - 0.9, y_gene[i], genes[i], adj = c(1, 0.5), cex = 0.64)
  for (i in seq_along(cell)) {
    text(x_dot[i], y_top + 1.5, cell[i], srt = 90, adj = c(0, 0.5),
         col = "black", cex = 0.57)
    text(109, mean(y_gene[((i-1)*5+1):(i*5)]), cell[i],
         adj = c(0, 0.5), col = lab$color[i], cex = 0.75, font = 2)
  }
  for (i in 1:25)
    segments(dot_l, y_top - i * 5 * pitch, dot_r, y_top - i * 5 * pitch,
             col = "#E1E1E1", lwd = 0.4)

  # 极浅的点阵网格与拟态表达点分层绘制；不把拟态当作作者测量值。
  grid_x <- rep(x_dot, times = length(genes))
  grid_y <- rep(y_gene, each = length(cell))
  points(grid_x, grid_y, pch = 16, cex = 0.055, col = "#DFE2E4")
  positive <- which(dots$pct_expressing > 0)
  col_index <- pmax(1, pmin(101, 1 + round(100 * dots$scaled_expression[positive])))
  points(x_dot[cell_index[positive]], y_gene[gene_index[positive]], pch = 16,
         cex = 0.11 + 0.54 * sqrt(dots$pct_expressing[positive]),
         col = cols[col_index])
  rect(dot_l, y_bottom, dot_r, y_top, border = "#363636", lwd = 0.7)

  # GO 条形仅占 0–3 的窄 NES 框；通路名称在框外，不挤进主图区域。
  bar_end <- nes_l + pmin(3, pmax(0, nes$nes)) / 3 * (nes_r - nes_l)
  rect(nes_l, y_go - pitch * 0.36, bar_end, y_go + pitch * 0.36,
       border = NA, col = lab$color[rep(seq_along(cell), each = 3)])
  rect(nes_l, y_bottom, nes_r, y_top, border = "#363636", lwd = 0.7)
  # 参考图的 NES 框右侧，逐条 GO 文字前各有一个短刻度。
  segments(nes_r, y_go, nes_r + 0.9, y_go, col = "#555555", lwd = 0.42)
  for (v in 0:3) {
    x <- nes_l + v / 3 * (nes_r - nes_l)
    segments(x, y_top, x, y_top + 1.5, lwd = 0.5)
    text(x, y_top + 3.2, v, cex = 0.72)
  }
  text((nes_l + nes_r)/2, y_top + 7, "NES", cex = 0.83)
  for (i in seq_len(nrow(nes)))
    text(go_x, y_go[i], nes$go_term[i], adj = c(0, 0.5), cex = 0.63)

  # 图例放在点图左上角的独立区域，色带紧靠上方标题。
  text(16, 344, "Expression", adj = c(0, 0.5), cex = 0.7, font = 3)
  text(16, 341, "scaled", adj = c(0, 0.5), cex = 0.7, font = 3)
  text(16, 338, "per-gene", adj = c(0, 0.5), cex = 0.7, font = 3)
  for (j in 1:101) {
    x0 <- 16 + (j-1)/101 * 25
    rect(x0, 334.0, x0 + 25/101, 336.3, border = NA, col = cols[j])
  }
  rect(16, 334.0, 41, 336.3, border = "grey45", lwd = 0.4)
  for (v in seq(0, 1, 0.2))
    text(16 + 25*v, 332.4, sprintf("%.1f", v), cex = 0.49)
  text(5, 333.4, "Percent", adj = c(0, 0.5), cex = 0.7, font = 3)
  text(5, 330.4, "of cells", adj = c(0, 0.5), cex = 0.7, font = 3)
  text(5, 327.4, "expressing", adj = c(0, 0.5), cex = 0.7, font = 3)
  for (i in 1:5) {
    pct <- rev(seq(0.2, 1, 0.2))[i]
    points(8.3, 324.2 - (i-1)*2.2, pch = 16, cex = 0.11+0.54*sqrt(pct), col = "grey55")
    text(11, 324.2 - (i-1)*2.2, sprintf("%d%%", round(pct*100)),
         adj = c(0, 0.5), cex = 0.53)
  }
  stopifnot(go_x + max(strwidth(nes$go_term, cex = 0.63)) < W)
}

# 四、PNG 和矢量 PDF 使用同一绘图函数，明确关闭设备并检查输出。
dir.create(file.path(root, "output"), showWarnings = FALSE)
png_path <- file.path(root, "output/preview.png")
pdf_path <- file.path(root, "output/preview.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "mm", res = 300, background = "white")
tryCatch(draw_panel(), finally = dev.off())
pdf(pdf_path, width = W/25.4, height = H/25.4, family = "Helvetica")
tryCatch(draw_panel(), finally = dev.off())
stopifnot(file.exists(png_path), file.exists(pdf_path),
          file.info(png_path)$size > 50000, file.info(pdf_path)$size > 10000)
cat("WROTE", png_path, file.info(png_path)$size, "bytes;", pdf_path,
    file.info(pdf_path)$size, "bytes\n")

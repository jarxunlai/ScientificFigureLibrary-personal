# Figure 167 / Figure 5 — Panel G 近似复刻
# cluster、作者标签、scMarkerAgent+LLM 标签来自公开补充表；
# positive/negative score 与 malignancy 为确定性虚拟构建，不是作者原始 UCell 数值。

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(grid)
})

# 一、输入与输出 ---------------------------------------------------------------
input_file <- "data/panelG_cluster_summary.csv"
out_dir <- "output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(out_dir, "panelG_recreated.png")
pdf_file <- file.path(out_dir, "panelG_recreated.pdf")

dat <- read_csv(input_file, show_col_types = FALSE) |>
  arrange(display_order)
stopifnot(nrow(dat) == 35L, identical(sort(as.integer(dat$cluster)), 0:34))

# 二、配色 --------------------------------------------------------------------
positive_col <- "#CB6C52"
negative_col <- "#364852"
malignancy_col <- c("Malignant cells" = "#541576", "Non-malignant cells" = "#DC461D")

model_col <- c(
  "astrocyte" = "#ADCBDE", "astrocyte-like cell" = "#AC9674", "B cell" = "#C5A2A4",
  "CD8-positive, alpha-beta cytotoxic T cell" = "#A490C0", "cDC2" = "#6E9A9F",
  "classical monocyte" = "#A4C88F", "cytotoxic T cell" = "#7B5FA2",
  "endothelial cell" = "#AC9674", "glioblastoma stem cell" = "#DA8686",
  "macrophage" = "#9CC671", "mesenchymal glioblastoma multiforme" = "#CD5152",
  "microglia/macrophage" = "#72AB50", "microglial cell" = "#729B4A",
  "neural progenitor cell" = "#7AA0C7", "neuron" = "#4F7DB3",
  "newly formed oligodendrocyte" = "#ADCBDE", "oligodendrocyte" = "#7AA0C7",
  "oligodendrocyte precursor-like cell" = "#4F7DB3", "pericyte" = "#D68B4C",
  "proneural" = "#C2292B", "radial glial cell" = "#ADCBDE",
  "tissue-resident memory T cells (T_rm)" = "#785E94",
  "tumor-associated macrophage" = "#BEB296", "Unknown" = "#CACACA",
  "vascular cell" = "#BEB296"
)

author_col <- c(
  "AC-like" = "#AC9674", "CD4/CD8" = "#C5A2A4", "DC" = "#6E9A9F",
  "Endothelial" = "#AC9674", "MES-like" = "#DA8686", "Mono" = "#A4C88F",
  "Mural cell" = "#D68B4C", "NPC-like" = "#CD5152",
  "Oligodendrocyte" = "#ADCBDE", "OPC-like" = "#C2292B", "RG" = "#7AA0C7",
  "TAM-BDM" = "#9CC671", "TAM-MG" = "#72AB50"
)

stopifnot(
  all(dat$model_label %in% names(model_col)),
  all(dat$author_label %in% names(author_col)),
  all(dat$malignancy %in% names(malignancy_col))
)

# 三、主图几何 -----------------------------------------------------------------
# 原图的 35 个 cluster 只占约 270°；左上四分之一是开口。
# 开口左端的五个竖向摘要槽与下方五条 annular track 共用同一组半径边界。
start_deg <- 90
end_deg <- -180
edges <- seq(start_deg, end_deg, length.out = nrow(dat) + 1)
centers <- (edges[-1] + edges[-length(edges)]) / 2

r_author <- c(18, 27)
r_model <- c(29, 38)
r_malig <- c(40, 49)
r_neg <- c(51, 65)
r_pos <- c(67, 81)

polar_xy <- function(cx, cy, r, deg) {
  rad <- deg * pi / 180
  c(cx + r * cos(rad), cy + r * sin(rad))
}

annular_polygon <- function(cx, cy, r0, r1, a0, a1, n = 10) {
  outer <- seq(a0, a1, length.out = n)
  inner <- seq(a1, a0, length.out = n)
  xo <- cx + r1 * cos(outer * pi / 180)
  yo <- cy + r1 * sin(outer * pi / 180)
  xi <- cx + r0 * cos(inner * pi / 180)
  yi <- cy + r0 * sin(inner * pi / 180)
  list(x = c(xo, xi), y = c(yo, yi))
}

draw_annular_cell <- function(cx, cy, r0, r1, a0, a1, fill = "white", col = "black", lwd = 0.55) {
  p <- annular_polygon(cx, cy, r0, r1, a0, a1)
  grid.polygon(unit(p$x, "mm"), unit(p$y, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd))
}

draw_radial_line <- function(cx, cy, r0, r1, deg, col, lwd = 0.8) {
  p0 <- polar_xy(cx, cy, r0, deg)
  p1 <- polar_xy(cx, cy, r1, deg)
  grid.lines(unit(c(p0[1], p1[1]), "mm"), unit(c(p0[2], p1[2]), "mm"), gp = gpar(col = col, lwd = lwd))
}

draw_main_open_ring <- function(cx, cy) {
  # 五条 annular tracks。
  for (i in seq_len(nrow(dat))) {
    a0 <- edges[i]
    a1 <- edges[i + 1]
    draw_annular_cell(cx, cy, r_author[1], r_author[2], a0, a1, author_col[[dat$author_label[i]]])
    draw_annular_cell(cx, cy, r_model[1], r_model[2], a0, a1, model_col[[dat$model_label[i]]])
    draw_annular_cell(cx, cy, r_malig[1], r_malig[2], a0, a1, malignancy_col[[dat$malignancy[i]]])
    draw_annular_cell(cx, cy, r_neg[1], r_neg[2], a0, a1, "white")
    draw_annular_cell(cx, cy, r_pos[1], r_pos[2], a0, a1, "white")

    # cluster id 位于 malignancy track。
    label_p <- polar_xy(cx, cy, mean(r_malig), centers[i])
    grid.text(
      as.character(dat$cluster[i]), unit(label_p[1], "mm"), unit(label_p[2], "mm"),
      rot = centers[i] - 270,
      gp = gpar(
        fontfamily = "Arial", fontsize = 5.2,
        col = if (dat$malignancy[i] == "Malignant cells") "white" else "black"
      )
    )

    # cluster-level score median [IQR]；扇区中间画 radial whisker + IQR 横杠。
    pos_base <- r_pos[1] + 0.8
    pos_q25 <- pos_base + dat$positive_q25[i] / 0.32 * (diff(r_pos) - 2.0)
    pos_med <- pos_base + dat$positive_median[i] / 0.32 * (diff(r_pos) - 2.0)
    pos_q75 <- pos_base + dat$positive_q75[i] / 0.32 * (diff(r_pos) - 2.0)
    draw_radial_line(cx, cy, pos_base, pos_q75, centers[i], positive_col, 0.55)
    for (rr in c(pos_q25, pos_med, pos_q75)) {
      p <- polar_xy(cx, cy, rr, centers[i])
      tang <- (centers[i] + 90) * pi / 180
      half <- if (rr == pos_med) 1.9 else 1.3
      grid.lines(
        unit(c(p[1] - half * cos(tang), p[1] + half * cos(tang)), "mm"),
        unit(c(p[2] - half * sin(tang), p[2] + half * sin(tang)), "mm"),
        gp = gpar(col = positive_col, lwd = if (rr == pos_med) 1.3 else 0.7)
      )
    }

    neg_base <- r_neg[1] + 0.8
    neg_q25 <- neg_base + dat$negative_q25[i] / 0.82 * (diff(r_neg) - 2.0)
    neg_med <- neg_base + dat$negative_median[i] / 0.82 * (diff(r_neg) - 2.0)
    neg_q75 <- neg_base + dat$negative_q75[i] / 0.82 * (diff(r_neg) - 2.0)
    draw_radial_line(cx, cy, neg_base, neg_q75, centers[i], negative_col, 0.55)
    for (rr in c(neg_q25, neg_med, neg_q75)) {
      p <- polar_xy(cx, cy, rr, centers[i])
      tang <- (centers[i] + 90) * pi / 180
      half <- if (rr == neg_med) 1.9 else 1.3
      grid.lines(
        unit(c(p[1] - half * cos(tang), p[1] + half * cos(tang)), "mm"),
        unit(c(p[2] - half * sin(tang), p[2] + half * sin(tang)), "mm"),
        gp = gpar(col = negative_col, lwd = if (rr == neg_med) 1.3 else 0.7)
      )
    }
  }

  # 圆环开口左端：五条 track 在同一 x 边界上向上延伸成竖向摘要槽。
  x_bounds <- c(
    cx - r_pos[2], cx - r_neg[2], cx - r_malig[2],
    cx - r_model[2], cx - r_author[2], cx - r_author[1]
  )
  y0 <- cy
  y1 <- cy + 61

  # 5 个摘要槽，左右边界正是开口处对应半径的端点。
  for (j in seq_len(5)) {
    grid.rect(
      x = unit(mean(x_bounds[j:(j + 1)]), "mm"), y = unit(mean(c(y0, y1)), "mm"),
      width = unit(abs(diff(x_bounds[j:(j + 1)])), "mm"), height = unit(y1 - y0, "mm"),
      gp = gpar(fill = "white", col = "black", lwd = 0.65)
    )
  }

  # Positive summary：使用虚拟分组值模拟原图箱线图方向。
  pos_mal <- dat$positive_median[dat$malignancy == "Malignant cells"]
  pos_non <- dat$positive_median[dat$malignancy == "Non-malignant cells"]
  neg_mal <- dat$negative_median[dat$malignancy == "Malignant cells"]
  neg_non <- dat$negative_median[dat$malignancy == "Non-malignant cells"]

  draw_box_summary <- function(xc, y_min, y_max, values, fill, width = 3.0) {
    q <- quantile(values, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
    scale_y <- function(v) y_min + v * (y_max - y_min)
    yq <- scale_y(q)
    grid.lines(unit(c(xc, xc), "mm"), unit(c(yq[1], yq[5]), "mm"), gp = gpar(col = "black", lwd = 0.55, lty = 2))
    grid.lines(unit(c(xc - width / 2, xc + width / 2), "mm"), unit(c(yq[1], yq[1]), "mm"), gp = gpar(col = "black", lwd = 0.55))
    grid.lines(unit(c(xc - width / 2, xc + width / 2), "mm"), unit(c(yq[5], yq[5]), "mm"), gp = gpar(col = "black", lwd = 0.55))
    grid.rect(unit(xc, "mm"), unit(mean(yq[2:4]), "mm"), unit(width, "mm"), unit(yq[4] - yq[2], "mm"), gp = gpar(fill = fill, col = "black", lwd = 0.65))
    grid.lines(unit(c(xc - width / 2, xc + width / 2), "mm"), unit(c(yq[3], yq[3]), "mm"), gp = gpar(col = "black", lwd = 0.75))
  }

  # Positive track summary 槽（最左）。
  pos_left <- x_bounds[1]
  pos_right <- x_bounds[2]
  draw_box_summary(pos_left + 0.34 * (pos_right - pos_left), cy + 7, cy + 44, pos_mal / 0.32, malignancy_col[["Malignant cells"]])
  draw_box_summary(pos_left + 0.72 * (pos_right - pos_left), cy + 7, cy + 44, pos_non / 0.32, malignancy_col[["Non-malignant cells"]])

  # Negative track summary 槽。
  neg_left <- x_bounds[2]
  neg_right <- x_bounds[3]
  draw_box_summary(neg_left + 0.34 * (neg_right - neg_left), cy + 7, cy + 44, neg_mal / 0.82, negative_col)
  draw_box_summary(neg_left + 0.72 * (neg_right - neg_left), cy + 7, cy + 44, neg_non / 0.82, malignancy_col[["Non-malignant cells"]])

  # Malignancy、模型标签、作者标签使用与原图一致的纵向 stacked strip。
  mal_counts <- prop.table(table(factor(dat$malignancy, levels = names(malignancy_col))))
  y <- cy
  for (name in names(malignancy_col)) {
    h <- as.numeric(mal_counts[[name]]) * (y1 - y0)
    grid.rect(unit(mean(x_bounds[3:4]), "mm"), unit(y + h / 2, "mm"),
              unit(abs(x_bounds[4] - x_bounds[3]), "mm"), unit(h, "mm"),
              gp = gpar(fill = malignancy_col[[name]], col = NA))
    y <- y + h
  }

  draw_stacked_strip <- function(left, right, labels, colors) {
    runs <- rle(labels)
    total <- length(labels)
    y <- cy
    idx <- 1
    for (k in seq_along(runs$lengths)) {
      h <- runs$lengths[k] / total * (y1 - y0)
      grid.rect(unit(mean(c(left, right)), "mm"), unit(y + h / 2, "mm"),
                unit(abs(right - left), "mm"), unit(h, "mm"),
                gp = gpar(fill = colors[[runs$values[k]]], col = "grey45", lwd = 0.25))
      y <- y + h
      idx <- idx + runs$lengths[k]
    }
  }
  draw_stacked_strip(x_bounds[4], x_bounds[5], dat$model_label, model_col)
  draw_stacked_strip(x_bounds[5], x_bounds[6], dat$author_label, author_col)

  # 摘要槽标签和显著性标记。
  labels <- c("Positive score", "Negative score", "Malignancy", "Cell type\n(scMarkerAgent+LLM)", "Cell type\n(author)")
  for (j in seq_along(labels)) {
    grid.text(
      labels[j], unit(mean(x_bounds[j:(j + 1)]), "mm"), unit(y1 + if (j <= 3) 8 else 5, "mm"),
      rot = 45, gp = gpar(fontfamily = "Arial", fontsize = if (j <= 3) 6.8 else 5.6)
    )
  }
  grid.text("***", unit(mean(x_bounds[1:2]), "mm"), unit(y1 - 4, "mm"), gp = gpar(fontfamily = "Arial", fontsize = 9))
  grid.text("***", unit(mean(x_bounds[2:3]), "mm"), unit(y1 - 4, "mm"), gp = gpar(fontfamily = "Arial", fontsize = 9))

  # score track 的刻度，仅在上端点附近显示，呼应原图。
  for (tick in c(0, 0.1, 0.2, 0.3)) {
    yy <- cy + tick / 0.32 * 14
    grid.text(format(tick, trim = TRUE), unit(cx - r_pos[2] - 2.2, "mm"), unit(yy, "mm"), just = "right", gp = gpar(fontfamily = "Arial", fontsize = 5.3))
  }
  for (tick in c(0, 0.2, 0.4, 0.6, 0.8)) {
    yy <- cy + tick / 0.82 * 14
    grid.text(format(tick, trim = TRUE), unit(cx - r_neg[2] - 1.8, "mm"), unit(yy, "mm"), just = "right", gp = gpar(fontfamily = "Arial", fontsize = 5.3))
  }
}

# 四、图例 --------------------------------------------------------------------
header_box <- function(label, fill, fontsize = 8.4) {
  grid.rect(gp = gpar(fill = fill, col = NA))
  grid.text(label, gp = gpar(fontfamily = "Arial", fontsize = fontsize))
}

item_grid <- function(labels, colors, ncol = 1, fontsize = 7.1) {
  n <- length(labels)
  nrow <- ceiling(n / ncol)
  pushViewport(viewport(layout = grid.layout(nrow, ncol)))
  for (i in seq_len(n)) {
    col_id <- floor((i - 1) / nrow) + 1
    row_id <- ((i - 1) %% nrow) + 1
    pushViewport(viewport(layout.pos.row = row_id, layout.pos.col = col_id))
    grid.rect(x = unit(1.8, "mm"), width = unit(2.5, "mm"), height = unit(2.5, "mm"),
              gp = gpar(fill = colors[[labels[[i]]]], col = "black", lwd = 0.55))
    grid.text(labels[[i]], x = unit(3.7, "mm"), just = "left",
              gp = gpar(fontfamily = "Arial", fontsize = fontsize))
    popViewport()
  }
  popViewport()
}

score_items <- function() {
  grid.circle(x = unit(2.0, "mm"), y = unit(0.68, "npc"), r = unit(1.7, "mm"), gp = gpar(fill = positive_col, col = NA))
  grid.text("Positive markers — median [IQR]", x = unit(4.7, "mm"), y = unit(0.68, "npc"), just = "left", gp = gpar(fontfamily = "Arial", fontsize = 7.5))
  grid.circle(x = unit(2.0, "mm"), y = unit(0.28, "npc"), r = unit(1.7, "mm"), gp = gpar(fill = negative_col, col = NA))
  grid.text("Negative markers — median [IQR]", x = unit(4.7, "mm"), y = unit(0.28, "npc"), just = "left", gp = gpar(fontfamily = "Arial", fontsize = 7.5))
}

# 五、总图 --------------------------------------------------------------------
draw_panel <- function() {
  grid.newpage()

  grid.text("G", unit(8, "mm"), unit(158, "mm"), just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = 15, fontface = "bold"))

  # 同一个主绘图坐标系包含开口圆环与左端五个竖向摘要槽。
  draw_main_open_ring(cx = 96, cy = 70)

  pushViewport(viewport(x = unit(205, "mm"), y = unit(149, "mm"), width = unit(68, "mm"), height = unit(7, "mm")))
  header_box("Malignancy signature score (UCell)", "#A20543", fontsize = 8.0)
  popViewport()
  pushViewport(viewport(x = unit(205, "mm"), y = unit(133, "mm"), width = unit(68, "mm"), height = unit(24, "mm")))
  score_items()
  popViewport()

  pushViewport(viewport(x = unit(279, "mm"), y = unit(149, "mm"), width = unit(54, "mm"), height = unit(7, "mm")))
  header_box("Malignancy", "#F36E43")
  popViewport()
  pushViewport(viewport(x = unit(279, "mm"), y = unit(133, "mm"), width = unit(54, "mm"), height = unit(24, "mm")))
  item_grid(names(malignancy_col), malignancy_col, ncol = 1, fontsize = 7.5)
  popViewport()

  pushViewport(viewport(x = unit(260, "mm"), y = unit(113, "mm"), width = unit(124, "mm"), height = unit(7, "mm")))
  header_box("Cell type (scMarkerAgent+LLM)", "#7FCBA3")
  popViewport()
  pushViewport(viewport(x = unit(260, "mm"), y = unit(81, "mm"), width = unit(124, "mm"), height = unit(58, "mm")))
  item_grid(names(model_col), model_col, ncol = 2, fontsize = 5.6)
  popViewport()

  pushViewport(viewport(x = unit(260, "mm"), y = unit(49, "mm"), width = unit(124, "mm"), height = unit(7, "mm")))
  header_box("Cell type (author)", "#FBDA82")
  popViewport()
  pushViewport(viewport(x = unit(260, "mm"), y = unit(32, "mm"), width = unit(124, "mm"), height = unit(27, "mm")))
  item_grid(names(author_col), author_col, ncol = 3, fontsize = 5.8)
  popViewport()
}

png(png_file, width = 340, height = 165, units = "mm", res = 300, type = "cairo", bg = "white")
draw_panel()
dev.off()

cairo_pdf(pdf_file, width = 340 / 25.4, height = 165 / 25.4, bg = "white")
draw_panel()
dev.off()

stopifnot(file.exists(png_file), file.info(png_file)$size > 0,
          file.exists(pdf_file), file.info(pdf_file)$size > 0)
message("Wrote ", png_file)
message("Wrote ", pdf_file)

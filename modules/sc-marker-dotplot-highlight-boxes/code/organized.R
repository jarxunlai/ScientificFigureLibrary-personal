# =============================================================================
# Nat Immunol Fig.1b 风格：胸腺 DC marker 气泡图 + 突出框
# =============================================================================
# 对照 Srinivasan et al., Nat Immunol 2026 Fig.1b（DOI 10.1038/s41590-025-02371-9）。
# 与上一版差异：x 轴用 11 个命名 cluster，不再用 Seurat 0–11；
# 三个框按文献列：cDC1 列 1–5 × Xcr1–Nr4a1；cDC2 列 7–9 × Sirpa–Cx3cr1；
# aDC 列 10–11 × Ccr7–Ccl22。数据仍为合成 DotPlot 长表。
# Community Archive 入口：Rscript render.R --input-dir <payload/data> --output <preview.png>
# 无 CLI 参数时仍读本条目 data/dotplot_dt.csv，并写 output/figures/。
# =============================================================================

library(ggplot2)
library(readr)
library(ragg)

parse_cli <- function(args) {
  input_dir <- NULL
  output_png <- NULL
  i <- 1L
  while (i <= length(args)) {
    if (identical(args[[i]], "--input-dir") && i < length(args)) {
      input_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(args[[i]], "--output") && i < length(args)) {
      output_png <- args[[i + 1L]]
      i <- i + 2L
    } else {
      stop("unexpected argument: ", args[[i]], call. = FALSE)
    }
  }
  list(input_dir = input_dir, output_png = output_png)
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
cli_mode <- !is.null(cli$input_dir) || !is.null(cli$output_png)
if (xor(is.null(cli$input_dir), is.null(cli$output_png))) {
  stop("usage: render.R --input-dir <dir> --output <preview.png>", call. = FALSE)
}

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
draft_dir <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir

if (cli_mode) {
  input_dir <- normalizePath(cli$input_dir, winslash = "/", mustWork = TRUE)
  csv_files <- list.files(input_dir, pattern = "\\.[Cc][Ss][Vv]$", full.names = TRUE)
  if (length(csv_files) != 1L) {
    stop("expected exactly one CSV in --input-dir", call. = FALSE)
  }
  input_csv <- csv_files[[1]]
  output_png <- cli$output_png
  output_pdf <- NULL
} else {
  input_csv <- file.path(draft_dir, "data", "dotplot_dt.csv")
  out_dir <- file.path(draft_dir, "output", "figures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  output_png <- file.path(out_dir, "marker_dotplot_boxes.png")
  output_pdf <- file.path(out_dir, "marker_dotplot_boxes.pdf")
}

plot_font_family <- "Arial"

dt <- as.data.frame(read_csv(input_csv, show_col_types = FALSE))
genes <- c(
  "Itgax", "Flt3", "Xcr1", "Irf8", "Cd36", "Cd207", "Nr4a1",
  "Ly6D", "Sirpa", "Irf4", "Tcf4", "Cd209a", "Mgl2", "Ccr2",
  "Epcam", "Ccl17", "Cd14", "Cx3cr1", "Ccr7", "Cd63", "Il12b",
  "Fabp5", "Ccl22"
)
clusters <- c(
  "Cycling (S)", "Cycling (G2M+S)", "Cycling (G2M) cDC1",
  "CD207hi Nur77-", "CD207lo Nur77+", "Ly6D+ cDC",
  "EpCAM+", "tDC2", "Cycling (G2M) cDC2", "aDC1", "aDC2"
)
dt$id <- factor(dt$id, levels = clusters)
dt$features.plot <- factor(dt$features.plot, levels = rev(genes))

p <- ggplot(dt, aes(x = id, y = features.plot)) +
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp), color = "black", shape = 21, stroke = 0.5) +
  xlab("") + ylab("") +
  scale_fill_gradientn(
    colours = c("#01009c", "#0000de", "#9559c8", "#faea4d", "#f09b37", "#ca2b1c", "#8b1a10"),
    limits = c(-2, 2), name = "Average expression"
  ) +
  scale_size(range = c(0, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = "Percent expressed")

p1 <- p +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9, color = "black", family = plot_font_family),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 11, face = "italic", color = "black", family = plot_font_family),
    legend.position = "right",
    legend.title = element_text(size = 11, family = plot_font_family)
  ) +
  guides(
    size = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
    fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

# y 因子从下到上是 Ccl22 → Itgax。文献三个框：
#   cDC1 列1–5 × Xcr1–Nr4a1；cDC2 列7–9 × Sirpa–Cx3cr1；aDC 列10–11 × Ccr7–Ccl22
p3 <- p1 +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
  annotate("rect", xmin = 0.5, xmax = 5.5, ymin = 16.5, ymax = 21.5, fill = NA, linewidth = 0.4, color = "black") +
  annotate("rect", xmin = 6.5, xmax = 9.5, ymin = 5.5, ymax = 15.5, fill = NA, linewidth = 0.4, color = "black") +
  annotate("rect", xmin = 9.5, xmax = 11.5, ymin = 0.5, ymax = 5.5, fill = NA, linewidth = 0.4, color = "black")

ggsave(output_png, p3, width = 8.2, height = 8.6, dpi = 300, bg = "white", device = ragg::agg_png)
if (!is.null(output_pdf)) {
  ggsave(output_pdf, p3, width = 8.2, height = 8.6, device = cairo_pdf, bg = "white")
}
message("wrote ", output_png)

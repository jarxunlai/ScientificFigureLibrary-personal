# =============================================================================
# Venn图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 整理相对原文：输出改到 output/figures；固定随机种子；关掉 VennDiagram 日志；
# 用填色三集合图作为 preview.png。作者 venn.diagram 参数原样保留。
# =============================================================================

# 优先：当前工作目录已是条目根（含 figure.yml）；否则用本脚本所在 code/。
# 避免 Rscript 包装器的 --file= 把 root 指到 _inventory。
if (file.exists(file.path(getwd(), "figure.yml")) &&
    file.exists(file.path(getwd(), "code", "organized.R"))) {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  this_file <- NULL
  frames <- sys.frames()
  if (length(frames)) {
    ofiles <- vapply(frames, function(e) {
      if (exists("ofile", envir = e, inherits = FALSE)) {
        as.character(get("ofile", envir = e, inherits = FALSE))
      } else NA_character_
    }, character(1))
    ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
    if (length(ofiles)) this_file <- ofiles[[length(ofiles)]]
  }
  if (is.null(this_file)) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) this_file <- sub("^--file=", "", file_arg)
  }
  script_dir <- if (!is.null(this_file)) dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)) else normalizePath(getwd())
  root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
}
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(grid)
  library(VennDiagram)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# 1. 运行设置
# 目的：可重复的模拟集合；日志写到输出目录后删除，避免污染条目根目录。
# -----------------------------------------------------------------------------
set.seed(1)
if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
}

# -----------------------------------------------------------------------------
# 2. 模拟三组集合
# 下游：3/4/5 维韦恩图共用这批抽样；检查点是各集合长度均为 200。
# -----------------------------------------------------------------------------
set1 <- paste(rep("word_", 200), sample(c(1:1000), 200, replace = FALSE), sep = "")
set2 <- paste(rep("word_", 200), sample(c(1:1000), 200, replace = FALSE), sep = "")
set3 <- paste(rep("word_", 200), sample(c(1:1000), 200, replace = FALSE), sep = "")

# -----------------------------------------------------------------------------
# 3. 基础三集合韦恩图
# -----------------------------------------------------------------------------
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1", "Set 2 ", "Set 3"),
  filename = file.path(out_dir, "venn_plot1.png"),
  output = TRUE,
  imagetype = "png",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw"
)

# -----------------------------------------------------------------------------
# 4. 填色、去描边的三集合（条目 preview）
# -----------------------------------------------------------------------------
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1", "Set 2 ", "Set 3"),
  filename = file.path(out_dir, "venn_plot2.png"),
  output = TRUE,
  imagetype = "png",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = "blank",
  fill = myCol,
  cex = 0.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  rotation = 1
)

# -----------------------------------------------------------------------------
# 5. 四集合、五集合
# -----------------------------------------------------------------------------
set4 <- paste(rep("word_", 200), sample(c(1:1000), 200, replace = FALSE), sep = "")
set5 <- paste(rep("word_", 200), sample(c(1:1000), 200, replace = FALSE), sep = "")

venn.diagram(
  x = list(set1, set2, set3, set4),
  category.names = c("Set 1", "Set 2 ", "Set 3", "Set 4"),
  filename = file.path(out_dir, "venn_plot3.png"),
  imagetype = "png",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw",
  col = "white",
  lty = 1,
  lwd = 1,
  fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0"),
  alpha = 0.90,
  label.col = "black",
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b"),
  cat.cex = 0.6,
  cat.fontfamily = "serif"
)

venn.diagram(
  x = list(set1, set2, set3, set4, set5),
  category.names = c("Set 1", "Set 2 ", "Set 3", "Set 4", "Set 5"),
  filename = file.path(out_dir, "venn_plot4.png"),
  imagetype = "png",
  height = 1000,
  width = 1000,
  resolution = 300,
  compression = "lzw",
  col = "white",
  lty = 1,
  lwd = 1,
  fill = c("#ffd7d8", "#d8f2e7", "#d9e7f2", "#eadff0", "#fff2cd"),
  alpha = 0.90,
  label.col = "black",
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#cb6274", "#7ba498", "#687d94", "#81668b", "#ffcf5c"),
  cat.cex = 0.6,
  cat.fontfamily = "serif",
  cat.pos = c(0, -30, -130, 130, 40),
  cat.dist = 0.18
)

# -----------------------------------------------------------------------------
# 6. 预览与清理
# 预览用填色三集合；删除 VennDiagram 写在输出目录的日志。
# -----------------------------------------------------------------------------
invisible(file.copy(
  file.path(out_dir, "venn_plot2.png"),
  file.path(root, "preview.png"),
  overwrite = TRUE
))
invisible(file.remove(list.files(out_dir, pattern = "\\.log$", full.names = TRUE)))
invisible(file.remove(list.files(root, pattern = "\\.log$", full.names = TRUE)))

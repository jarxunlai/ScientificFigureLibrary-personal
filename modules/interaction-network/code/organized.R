# =============================================================================
# 相互作用网络图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 整理相对原文：每张 igraph 图单独 png 保存；HTML widget 写到条目 HtmlWidget/；
# 用 mtcars 高相关网络作为 preview.png。作者 plot() 参数原样保留。
# =============================================================================

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
html_dir <- file.path(root, "HtmlWidget")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(igraph)
  library(RColorBrewer)
  library(networkD3)
  library(htmlwidgets)
})

# -----------------------------------------------------------------------------
# 1. 四种 igraph 布局对照
# 目的：同一随机邻接矩阵上比较 sphere / circle / random / FR。
# 检查点：20 个节点；边由 0.2 的抽样概率产生，布局只影响位置。
# -----------------------------------------------------------------------------
set.seed(1)
data <- matrix(sample(0:1, 400, replace = TRUE, prob = c(0.8, 0.2)), nrow = 20)
data <- pmax(data, t(data))
diag(data) <- 0
network <- graph_from_adjacency_matrix(data, mode = "undirected", diag = FALSE)

png(file.path(out_dir, "igraph_layouts.png"), width = 1600, height = 1600, res = 150, bg = "white")
opar <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1), bg = "white")
plot(network, layout = layout.sphere, main = "sphere")
plot(network, layout = layout.circle, main = "circle")
plot(network, layout = layout.random, main = "random")
plot(network, layout = layout.fruchterman.reingold, main = "fruchterman.reingold")
par(opar)
dev.off()

# -----------------------------------------------------------------------------
# 2. 修改顶点 / 边外观
# -----------------------------------------------------------------------------
set.seed(1)
data <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.8, 0.2)), nc = 10)
data <- pmax(data, t(data))
diag(data) <- 0
network <- graph_from_adjacency_matrix(data, mode = "undirected", diag = FALSE)

png(file.path(out_dir, "igraph_styled.png"), width = 1400, height = 1400, res = 150, bg = "black")
par(bg = "black", mar = c(1, 1, 1, 1))
plot(
  network,
  vertex.color = rgb(0.8, 0.4, 0.3, 0.8),
  vertex.frame.color = "white",
  vertex.shape = "circle",
  vertex.size = 14,
  vertex.size2 = NA,
  vertex.label = LETTERS[1:10],
  vertex.label.color = "white",
  vertex.label.family = "serif",
  vertex.label.font = 2,
  vertex.label.cex = 1,
  vertex.label.dist = 0,
  vertex.label.degree = 0,
  edge.color = "white",
  edge.width = 4,
  edge.arrow.size = 1,
  edge.arrow.width = 1,
  edge.lty = "solid",
  edge.curved = 0.3
)
dev.off()

# -----------------------------------------------------------------------------
# 3. mtcars 高相关网络（条目 preview）
# 节点颜色映射气缸数 cyl。
# -----------------------------------------------------------------------------
mat <- cor(t(mtcars[, c(1, 3:6)]))
mat[mat < 0.995] <- 0
mat <- pmax(mat, t(mat))
diag(mat) <- 0
network <- graph_from_adjacency_matrix(mat, weighted = TRUE, mode = "undirected", diag = FALSE)
coul <- brewer.pal(nlevels(as.factor(mtcars$cyl)), "Set2")
my_color <- coul[as.numeric(as.factor(mtcars$cyl))]

png(file.path(out_dir, "mtcars_network.png"), width = 1400, height = 1400, res = 150, bg = "grey13")
par(bg = "grey13", mar = c(0, 0, 0, 0))
set.seed(4)
plot(
  network,
  vertex.size = 12,
  vertex.color = my_color,
  vertex.label.cex = 0.7,
  vertex.label.color = "white",
  vertex.frame.color = "transparent"
)
text(0, 0, "mtcars network", col = "white", cex = 1.5)
legend(
  x = -0.2, y = -0.12,
  legend = paste(levels(as.factor(mtcars$cyl)), " cylinders", sep = ""),
  col = coul,
  bty = "n", pch = 20, pt.cex = 2, cex = 1,
  text.col = "white", horiz = FALSE
)
dev.off()

# -----------------------------------------------------------------------------
# 4. networkD3 交互式网络
# HTML 不是 PNG preview；selfcontained=FALSE，避免依赖 pandoc。
# -----------------------------------------------------------------------------
data <- data.frame(
  from = c("A", "A", "B", "D", "C", "D", "E", "B", "C", "D", "K", "A", "M"),
  to = c("B", "E", "F", "A", "C", "A", "B", "Z", "A", "C", "A", "B", "K"),
  stringsAsFactors = FALSE
)

p <- simpleNetwork(data, height = "400px", width = "400px")
saveWidget(
  p,
  file = file.path(html_dir, "networkInteractive1.html"),
  selfcontained = FALSE
)

p <- simpleNetwork(
  data,
  height = "400px",
  width = "400px",
  Source = 1,
  Target = 2,
  linkDistance = 10,
  charge = -900,
  fontSize = 14,
  fontFamily = "serif",
  linkColour = "#666",
  nodeColour = "#69b3a2",
  opacity = 0.9,
  zoom = TRUE
)
saveWidget(
  p,
  file = file.path(html_dir, "networkInteractive2.html"),
  selfcontained = FALSE
)

# -----------------------------------------------------------------------------
# 5. 预览：带气缸着色与图例的相关网络
# -----------------------------------------------------------------------------
invisible(file.copy(
  file.path(out_dir, "mtcars_network.png"),
  file.path(root, "preview.png"),
  overwrite = TRUE
))
while (grDevices::dev.cur() > 1) grDevices::dev.off()

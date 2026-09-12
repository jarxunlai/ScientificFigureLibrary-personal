# =============================================================================
# 层状结构可视化
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 整理相对原文：三张 ggraph 都 ggsave；旭日图 CSV 走 data/；环形树作 preview。
# 作者布局与 aes 原样保留，不改成自写 helper。
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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(ggraph)
  library(igraph)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# 1. 模拟三级层次，画环形树
# 目的：origin → group0* → group* → subgroup_*；叶子按组排序后绕圆排布。
# 检查点：set.seed(1) 后边表可复现；这张作为条目 preview。
# -----------------------------------------------------------------------------
set.seed(1)
d0 <- data.frame(from = "origin", to = paste("group0", seq(1, 3), sep = ""))
d1 <- data.frame(from = rep(d0$to, c(5, 3, 2)), to = paste("group", seq(1, 10), sep = ""))
d2 <- data.frame(
  from = sort(sample(rep(d1$to, each = 100), 100, replace = FALSE)),
  to = sample(paste("subgroup", seq(1, 100), sep = "_"), 100, replace = FALSE)
)

group_sum <- summary(as.factor(d2$from))
d1$group_sum <- group_sum[d1$to]
d1 <- d1 %>%
  mutate(order1 = as.numeric(factor(
    from,
    levels = unique(from)[sort(summary(as.factor(from)), index.return = TRUE, decreasing = TRUE)$ix],
    order = TRUE
  ))) %>%
  arrange(order1, desc(group_sum))
d2 <- left_join(d2, d1, by = c("from" = "to"))
d2 <- d2 %>%
  mutate(order2 = as.numeric(factor(
    from,
    levels = unique(from)[sort(summary(as.factor(from)), index.return = TRUE, decreasing = TRUE)$ix],
    order = TRUE
  ))) %>%
  arrange(order1, desc(group_sum), order2)
edges <- rbind(d0, d1[, 1:2], d2[, 1:2])

vertices_name <- unique(c(as.character(edges$from), as.character(edges$to)))
vertices <- data.frame(name = vertices_name, value = runif(length(vertices_name)))
rownames(vertices) <- vertices_name

d2 <- d2 %>%
  left_join(vertices, by = c("to" = "name")) %>%
  arrange(order1, desc(group_sum), order2, desc(value))
edges <- rbind(d0, d1[, 1:2], d2[, 1:2])

list_unique <- unique(c(as.character(edges$from), as.character(edges$to)))
vertices <- data.frame(
  name = list_unique,
  value = vertices[list_unique, "value"]
)

vertices$group <- edges$from[match(vertices$name, edges$to)]
vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from)))
nleaves <- length(myleaves)
vertices$id <- nleaves / 360 * 90
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)

mygraph <- graph_from_data_frame(edges, vertices = vertices, directed = TRUE)

p_circular <- ggraph(mygraph, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(aes(colour = after_stat(index))) +
  scale_edge_colour_distiller(palette = "RdPu") +
  geom_node_text(aes(x = x * 1.25, y = y * 1.25, angle = angle, label = name, color = group),
                 size = 2.7, alpha = 1) +
  geom_node_point(aes(x = x * 1.07, y = y * 1.07, fill = group, size = value),
                  shape = 21, stroke = 0.2, color = "black", alpha = 1) +
  scale_colour_manual(values = rep(brewer.pal(9, "Paired"), 30)) +
  scale_fill_manual(values = rep(brewer.pal(9, "Paired"), 30)) +
  scale_size_continuous(range = c(0.1, 7)) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

ggsave(
  file.path(out_dir, "circular_dendrogram.png"),
  plot = p_circular,
  width = 8, height = 8, dpi = 200, bg = "white"
)

# -----------------------------------------------------------------------------
# 2. 旭日图：扇区等分
# 数据：data/旭日图.csv，季→月→周；外圈按 Value 填色，角度均分。
# -----------------------------------------------------------------------------
df <- read.csv(file.path(root, "data", "旭日图.csv"), header = TRUE, stringsAsFactors = FALSE)

edges <- data.frame(rbind(
  cbind(rep("origin", 4), unique(as.character(df$Season))),
  as.matrix(df[!duplicated(df[c("Season", "Month")]), 1:2]),
  as.matrix(df[!duplicated(df[c("Month", "Week")]), 2:3])
))
colnames(edges) <- c("from", "to")

vertices0 <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
df_leaf <- df[, c("Week", "Value", "label")]
df_leaf$angle <- 90 - (1:nrow(df_leaf)) / nrow(df_leaf) * 360
df_leaf$angle <- ifelse(df_leaf$angle < -90, df_leaf$angle + 180, df_leaf$angle)
vertices <- left_join(vertices0, df_leaf, by = c("name" = "Week"))
graph <- graph_from_data_frame(edges, vertices = vertices)

p_sun_equal <- ggraph(graph, layout = "partition", circular = TRUE) +
  geom_node_arc_bar(aes(filter = (depth <= 2 & depth > 0)), fill = "#FEE5D9", color = "black", linewidth = 0.1) +
  geom_node_arc_bar(aes(filter = (depth == 3), fill = Value), linewidth = 0.1) +
  geom_node_text(aes(filter = (depth <= 2 & depth > 0), label = name, size = -depth), angle = 0, colour = "black") +
  geom_node_text(aes(filter = !is.na(label), label = label, angle = angle), size = 2, colour = "black") +
  scale_size(range = c(3, 5)) +
  coord_fixed() +
  scale_fill_distiller(palette = "Reds") +
  guides(size = "none") +
  theme_void()

ggsave(
  file.path(out_dir, "sunburst_equal.png"),
  plot = p_sun_equal,
  width = 8, height = 8, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# 3. 旭日图：扇区与 Value 成比例
# 用重复叶子节点把 partition 角度拉成与 Value 近似成比例。
# -----------------------------------------------------------------------------
fake_circle <- c()
for (i in 1:nrow(df)) {
  fake_circle <- append(fake_circle, rep(df$Week[i], round(10 * df$Value[i])))
}

edges <- data.frame(rbind(
  cbind(rep("origin", 4), unique(as.character(df$Season))),
  as.matrix(df[!duplicated(df[c("Season", "Month")]), 1:2]),
  as.matrix(df[!duplicated(df[c("Month", "Week")]), 2:3]),
  cbind(fake_circle, as.character(1:length(fake_circle)))
))
colnames(edges) <- c("from", "to")

vertices0 <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))))
df_leaf <- df[, c("Week", "Value", "label")]
df_leaf$angle <- 90 - (cumsum(df_leaf$Value) - df_leaf$Value / 2) / sum(df_leaf$Value) * 360
df_leaf$angle <- ifelse(df_leaf$angle < -90, df_leaf$angle + 180, df_leaf$angle)
vertices <- left_join(vertices0, df_leaf, by = c("name" = "Week"))

df_color <- data.frame(rbind(
  as.matrix(df[!duplicated(df[c("Season", "Season")]), c(1, 1)]),
  as.matrix(df[!duplicated(df[c("Season", "Month")]), c(1, 2)]),
  as.matrix(df[!duplicated(df[c("Season", "Week")]), c(1, 3)])
))
colnames(df_color) <- c("Season", "name")
vertices <- left_join(vertices, df_color, by = "name")
graph <- graph_from_data_frame(edges, vertices = vertices)

p_sun_prop <- ggraph(graph, layout = "partition", circular = TRUE) +
  geom_node_arc_bar(aes(filter = (depth <= 3 & depth > 0), fill = Season), linewidth = 0.1) +
  geom_node_text(aes(filter = (depth <= 2 & depth > 0), label = name, size = -depth), angle = 0, colour = "black") +
  geom_node_text(aes(filter = !is.na(label), label = label, angle = angle), size = 2, colour = "black") +
  scale_size(range = c(3, 4)) +
  coord_fixed() +
  scale_fill_brewer(palette = "Reds", direction = -1) +
  guides(size = "none", fill = "none") +
  theme_void()

ggsave(
  file.path(out_dir, "sunburst_proportional.png"),
  plot = p_sun_prop,
  width = 8, height = 8, dpi = 300, bg = "white"
)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
invisible(file.copy(
  file.path(out_dir, "circular_dendrogram.png"),
  file.path(root, "preview.png"),
  overwrite = TRUE
))

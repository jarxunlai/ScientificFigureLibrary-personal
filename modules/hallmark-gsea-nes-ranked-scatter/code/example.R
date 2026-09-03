# 绘图层入口：读取通路级 GSEA 表并按作者脚本的 ggplot 映射作图。
# 不读取论文补充表，不调用 clusterProfiler::GSEA。
# 作者原脚本见 original.R。

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") {
  dirname(script_dir)
} else {
  script_dir
}

data <- read.csv(file.path(root, "data/example_gsea.csv"), stringsAsFactors = FALSE)
data$setSize_1 <- data$setSize / 10
data <- data[order(data$NES, decreasing = TRUE), ]
data$ID <- factor(data$ID, levels = data$ID)
data$xlab <- seq_len(nrow(data))
if (nrow(data) != 49) {
  warning("作者原脚本把 xlab 写成 1:49；当前表有 ", nrow(data), " 行。")
}

label <- c(
  "KRAS_SIGNALING_DN", "INTERFERON_ALPHA_RESPONSE",
  "UV_RESPONSE_DN", "EMT", "TNFA_SIGNALING_VIA_NFKB", "MYC_TARGETS_V2",
  "KRAS_SIGNALING_UP", "MYC_TARGETS_V1", "G2M_CHECKPOINT", "E2F_TARGETS"
)
data_label <- data[data$ID %in% label, ]
data_label$col <- c("black", "grey80", "grey80", "grey80", "grey80", "#ffb882", "grey80", "#ff5eff", "#ff5eff")

p <- ggplot(data = data, aes(x = xlab, y = NES)) +
  geom_point(
    aes(size = setSize_1, alpha = -log10(pvalue)),
    shape = 21, stroke = 0.7, fill = "#0000ff", colour = "black"
  ) +
  scale_size_continuous(range = c(0.2, 8)) +
  xlab(label = "Hallmark gene sets") +
  ylab(label = "Normalized enrichment score (NES)") +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = seq(0, 50, by = 10), labels = seq(0, 50, by = 10)) +
  scale_y_continuous(breaks = seq(-4, 2.3, by = 1), labels = seq(-4, 2.3, by = 1)) +
  guides(
    size = guide_legend(title = "Detection\n(proportion)"),
    alpha = guide_legend(title = "Significance\n(-log10 p-val.)")
  ) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 13)
  )

p3 <- p +
  geom_text_repel(
    data = data_label,
    aes(x = xlab, y = NES, label = ID),
    size = 3,
    color = data_label$col,
    force = 20,
    point.padding = 0.5,
    min.segment.length = 0,
    hjust = 1.2,
    segment.color = "grey20",
    segment.size = 0.3,
    segment.alpha = 0.8,
    nudge_y = -0.1
  )

out_png <- file.path(root, "preview.png")
ggsave(filename = out_png, plot = p3, width = 6.2, height = 5, dpi = 300, bg = "white")
ggsave(filename = "Figure1C.pdf", plot = p3, width = 6.2, height = 5, bg = "white")

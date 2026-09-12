# =============================================================================
# 形状火山图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# =============================================================================

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
set.seed(20260912)
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
})


library(latex2exp)
library(ggrepel)

# load the data
# 独立生成的模拟差异表；不代表真实分析结果。
deg_data <- data.frame(log2FoldChange = rnorm(1000, 0, 2), padj = 10^(-runif(1000, 0, 12)), row.names = paste0("Gene_", seq_len(1000)))

# 任意指定一组分组变量 -- 无实际意义
deg_data$group <- ifelse(sample(1:nrow(deg_data), replace = T) > nrow(deg_data)/2,
                         "Group A", "Group B")

# -log10(padj)小于2的设置为灰色
deg_data$color <- ifelse(-log10(deg_data$padj) < 2, "#bfc0c1",
                         ifelse(deg_data$log2FoldChange>0, "#e88182", "#6489b2") )

head(deg_data)
# row   baseMean log2FoldChange     lfcSE       stat       pvalue
# padj               group   color     label

# 绘图
ggplot(deg_data)+
  geom_point(aes(log2FoldChange, -log10(padj),
                 shape = group, size = -log10(padj)),
             color = deg_data$color,
             alpha = 0.7
             )+
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_hline(yintercept = 2, linetype = "longdash")+
  scale_size_continuous(range = c(0.2, 3))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = TeX("$Log_2 \\textit{FC}$"),
       y = expression(-log[10](FDR)))

primary_plot <- ggplot2::last_plot()
ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 6)


# 添加标签：
deg_data$label <- rep("", nrow(deg_data))
deg_data$label[order(deg_data$padj)[1:20]] <- rownames(deg_data)[order(deg_data$padj)[1:20]]


ggplot(deg_data, aes(log2FoldChange, -log10(padj)))+
  geom_point(aes(shape = group, size = -log10(padj)),
                 color = deg_data$color,
                 alpha = 0.7)+
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_hline(yintercept = 2, linetype = "longdash")+
  geom_text_repel(aes(label = label), size = 2, color = deg_data$color,
                  max.overlaps = 100, key_glyph = draw_key_point)+
  scale_size_continuous(range = c(0.2, 3))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x = TeX("$Log_2 \\textit{FC}$"),
       y = expression(-log[10](FDR)))

ggsave(file.path(out_dir, "plot2.pdf"), height = 5, width = 6)

ggplot2::ggsave(file.path(root, "preview.png"), plot = primary_plot, width = 6, height = 5, dpi = 200, bg = "white")

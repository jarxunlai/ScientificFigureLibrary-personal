# 固定模拟示例的随机种子，便于重复生成预览。
set.seed(20260911)
# =============================================================================
# 箱线图加热图注释
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

# library(tidyverse) 已在文件头拆包

# 构造数据：
# size_fractions:
Size_fractions <- as.data.frame(matrix(NA, nrow = 30, ncol = 5))
for (i in 1:5) {
  Size_fractions[,i] <- c(runif(10, min = 2, max = 60),
                         runif(10, min = 1, max = 30),
                         runif(10, min = 0, max = 5))
}
colnames(Size_fractions) <- c("<0.2","0.2-0.8","0.2-3","0.8-20",">0.2")

# Latitudes
Latitudes <- as.data.frame(matrix(NA, nrow = 30, ncol = 3))

for (i in 1:3) {
  Latitudes[,i] <- c(runif(10, min = 40, max = 80),
                     runif(10, min = 5, max = 60),
                     runif(10, min = 0, max = 70))
}

colnames(Latitudes) <- c("<30", "30-60", ">60")

# Depth layers
Depth_layers <- as.data.frame(matrix(NA, nrow = 30, ncol = 3))

for (i in 1:3) {
  Depth_layers[,i] <- c(runif(10, min = 40, max = 80),
                     runif(10, min = 5, max = 40),
                     runif(10, min = 0, max = 60))
}

colnames(Depth_layers) <- c("EPI", "MES", "BAT")

# 合并：
data <- cbind(Size_fractions, Latitudes, Depth_layers)
data$group <- factor(rep(c("OMD", "GEM", "GORG"), each = 10),
                     levels = c("OMD", "GEM", "GORG"))

# 转长数据：
data_long <- pivot_longer(data, cols = !group, 
                          names_to = "x", values_to = "y")
data_long$group2 <- factor(rep(rep(c("Size fraction(um)", "Latitudes", "Depth layers"),
                        c(5, 3, 3)), 30), levels = c("Size fraction(um)", "Latitudes", "Depth layers"))
data_long$x <- factor(data_long$x, levels = unique(data_long$x))


# 画图：
p1 <- ggplot(data_long)+
  geom_boxplot(aes(x, y, fill = group, color = group), alpha = 0.5, outlier.shape = NA)+
  facet_grid(.~group2, scales="free_x", space = "free_x")+
  scale_color_manual(name = "", values = c("#5aadd0", "#a9dcb6", "#af7aa1"))+
  scale_fill_manual(name = "", values = c("#5aadd0", "#a9dcb6", "#af7aa1"))+
  xlab("")+
  ylab("Reads mapped(%)")+
  scale_y_continuous(breaks = seq(0,80,20), limits = c(0, 80),expand = c(0,0))+
  # 主题调整：
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

p1

ggsave(file.path(out_dir, "plot1.pdf"), height = 3.5, width = 8)

# 对齐上方类别注释和下方分面箱线图；不改变数据、顺序或配色。
rect_data <- dplyr::distinct(data_long, x, group2)
p2 <- ggplot(rect_data, aes(x, 1)) +
  geom_tile(aes(fill = x), color = "white", linewidth = 1, width = 1, height = 1) +
  geom_text(aes(label = x), angle = 90, color = "white", size = 3) +
  facet_grid(. ~ group2, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#fa7f5e", "#f2605c", "#dd4968", "#c43c75", "#a8327d", "#d38c39", "#cac5bf", "#4f69b1", "#dbdb5f", "#509f9b", "#1e497f")) +
  theme_void() + theme(legend.position = "none", strip.text = element_blank(), panel.spacing.x = unit(5.5, "pt"))
library(patchwork)
p <- p2 / p1 + plot_layout(heights = c(1.2, 5))
ggsave(file.path(out_dir, "plots.pdf"), plot = p, height = 5, width = 8)

# 显式生成发布预览，不依赖外部图片转换。
ggplot2::ggsave(file.path(root, "preview.png"), plot = p, device = ragg::agg_png, width = 8, height = 5, units = "in", dpi = 300, bg = "white")

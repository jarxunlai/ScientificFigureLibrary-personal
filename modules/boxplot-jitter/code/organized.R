# 固定模拟示例的随机种子，便于重复生成预览。
set.seed(20260911)
# =============================================================================
# 箱线抖动散点图
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

########### 数据构建 --------------
G1 <- runif(100, min = 0, max = 7)
G2 <- runif(20, min = 5, max = 7)
G3 <- runif(10, min = 1, max = 6)
G4 <- runif(15, min = 2, max = 6)
G5 <- runif(20, min = 2.2, max = 6.5)
G6 <- runif(10, min = 3.5, max = 5)
G7 <- runif(80, min = 1, max = 6)
G8 <- runif(70, min = 1, max = 5.5)
G9 <- runif(60, min = 1.5, max = 6)
G10 <- runif(200, min = 1, max = 7.2)

# 合并：
data <- data.frame(Group = rep(paste0("G", 1:10),
                                c(100, 20, 10, 15,
                                  20, 10, 80, 70,
                                  60, 200)),
                   values = c(G1,G2,G3,G4,G5,G6,G7,G8,G9,G10))

data$Group <- factor(data$Group, levels = paste0("G", 1:10))

head(data)
#   Group     values
# 1    G1 2.20054387
# 2    G1 1.90207512
# 3    G1 2.74224843
# 4    G1 2.17059052
# 5    G1 4.14728737
# 6    G1 0.01258516


############# 绘图 -----------
library(ggplot2)
library(latex2exp)

ggplot(data, aes(Group, values))+
  # 箱线图：
  geom_boxplot(outlier.shape = NA, width = 0.6)+
  # 抖动散点：
  geom_jitter(aes(color = Group), width = 0.15, size = 1)+
  # 横线：
  geom_hline(yintercept = 4, linetype = "dashed")+
  # 箭头：
  geom_segment(aes(x = 2, y = 7.5, xend = 2, yend = 7.2),
               arrow = arrow(length = unit(1, "mm"))) +
  scale_color_manual(name = "Subtype",
                     values = c("#fd6ab0", "#aa5700", "#f48326", "#ffd711",
                                "#9bd53f", "#00ae4c", "#00c1e3", "#007ddb",
                                "#8538d1", "#d01910"))+
  # 文字注释：
  annotate("text", label = "P = 1.6e-06 (illustrative)",
           size = 3, x = 2, y = 8)+
  # 坐标轴标签：
  xlab("")+
  ylab(TeX("$Log_{2}(FPKM+1)$"))+
  # 标题：
  ggtitle(TeX("$\\textit{GATA3}$ gene expression in T-ALL"))+
  # 主题：
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(override.aes = list(size=2),
                            title.theme = element_text(face = "bold")))

ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 7)

# 显式生成发布预览，不依赖外部图片转换。
ggplot2::ggsave(file.path(root, "preview.png"), plot = last_plot(), device = ragg::agg_png, width = 7, height = 5, units = "in", dpi = 300, bg = "white")

# =============================================================================
# 散点等高线图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」33散点图+等高线
# 本地复现：drafts/ks-shengxin-huitu-repro/c33b-scatter-contour
# Pixi：项目根 default 环境；library(tidyverse) 已拆成 ggplot2/dplyr/tidyr 等。
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

# 载入R包：
library(ggplot2)
# library(tidyverse) 已在文件头拆包

# 创建数据
cluster1 <- data.frame(x = rnorm(5000, 6, 1.5),
                       y = rnorm(5000, 6, 1.4) )
cluster2 <- data.frame(x = rnorm(5000, 11, 1.6),
                       y = rnorm(5000, 7, 1.7) )
data <- rbind(cluster1, cluster2)
data$group <- rep(c("A", "B"), each = 5000)

# 绘制基础散点图
ggplot(data, aes(x=x, y=y)) +
  geom_point()

# 很显然，你不能明显的看出这是两坨不一样的散点，
# 但是加上等高线就不一样了
ggplot(data, aes(x=x, y=y)) +
  geom_point()+
  geom_density2d()

# 很明显这样就很明显可以看出两坨了
# 我们再来做一波美化:
ggplot(data, aes(x=x, y=y)) +
  # 绘制最底层散点图：
  geom_point(aes(color = group), size = 1)+
  # 绘制等高线的填充：
  geom_density_2d_filled(alpha = 0.3) +
  # 绘制等高线：
  geom_density2d(color = "#d14524") +
  # 调整散点颜色：
  scale_color_brewer(palette="Set2", direction=1)+
  # 调整填充色渐变：
  scale_fill_brewer(palette=4, direction=1)+
  # 设置x和y轴拓展：
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Feature1")+
  ylab("Feature2")+
  # 设置主题：
  theme_bw()+
  # 设置图例：
  guides(color = "none")

ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 7)

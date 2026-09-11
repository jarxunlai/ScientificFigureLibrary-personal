# =============================================================================
# 分组重叠直方图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」32分组重叠直方图
# 本地复现：drafts/ks-shengxin-huitu-repro/c32-overlap-histogram
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

########## 分组直方图 ---------
library(ggplot2)

# 构造数据：
data <- data.frame(Control = rnorm(1000, mean = 500, sd = 200),
                   Case = rnorm(1000, mean = 400, sd = 200))

head(data)
#    Control       Case
# 1 388.3990  -7.477409
# 2 504.1061 379.507800
# 3 448.9444 309.780209
# 4 434.8209 597.555137
# 5 757.1036 308.031864
# 6 741.7811 593.934408


####### 绘图 -------------
p <- ggplot(data)+
  geom_histogram(aes(Case), binwidth = 20,
                 color = "white", fill = "#cc7833", alpha = 0.5)+
  geom_histogram(aes(Control), binwidth = 20,
                 color = "white", fill = "#75aadb", alpha = 0.5)+
  ggtitle("Density Histogram")+
  xlab("Density")+
  ylab("Expression")+
  theme_bw()

p

# 但是由于没有分组变量，图例是不会自动生成的，怎么办？
# 创建一个空的数据框
tmp <- data.frame(x = c(0, 0), y = c(0, 0), group = c("A", "B"))

# 创建绘图对象并设置主题
p1 <- p + geom_bar(data = tmp, aes(x = x, y = y, fill = group),
                   stat = "identity")+
  scale_fill_manual(name = "Group", labels = c("Control", "Case"),
                    values = c("#cc7833", "#75aadb"))+
  scale_x_continuous(breaks = seq(-200, 1000,200))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6, hjust = 0.5),
        legend.box.background = element_rect(color = "#aaaaaa"))
p1

ggsave(file.path(out_dir, "plot.pdf"), height = 4, width = 5)

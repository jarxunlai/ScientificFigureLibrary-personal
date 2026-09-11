# =============================================================================
# 分半箱线抖动散点图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」041分半箱线图+抖动散点图
# 本地复现：drafts/ks-shengxin-huitu-repro/041-half-boxplot-jitter
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

# library(tidyverse) 已在文件头拆包
# devtools::install_local("gghalves-master.zip")
library(gghalves)

min <- runif(20, min = -0.3, max = 0.1)
max <- runif(20, min = 0.2, max = 0.9)

####### 构造模拟数据 ---------------
data <- data.frame(otus = rep(paste0("otu", 1), 100),
                   correlation = runif(100, min = min[1], max = max[1]),
                   group = sample(c("Specialist phage", "Generalist phage"),
                                  100, replace = T))

# 批量构建数据，这一步的目的是尽量保证组与组之间差别大一些
for (i in 2:20) {
  data_tmp <- data.frame(otus = rep(paste0("otu", i), 100),
                         correlation = runif(100, min = min[i], max = max[i]),
                         group = sample(c("Specialist phage", "Generalist phage"),
                                        100, replace = T))
  data <- rbind(data, data_tmp)
}

########### 绘图 -------------------
ggplot(data, aes(otus, correlation))+
  # 绘制分半箱线图：
  geom_half_boxplot(fill = "#75aadb", alpha = 0.6)+
  # 绘制分半散点图：
  geom_half_point_panel(aes(color = group), size = 0.5,
                        # 调整抖动散点的宽度：
                        position = position_jitter(width = 0.2))+
  # 添加横线：
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
  # 散点颜色模式：
  scale_color_manual(name = "", values = c("#a1b4c3", "#fbc47e"))+
  # 主题调整：
  theme_classic()+
  # 调整图例位置：
  theme(legend.position = "top")+
  # 修改图例散点大小：
  guides(color = guide_legend(override.aes = list(size = 4)))

# 保存：
ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 10)

p041 <- ggplot(data, aes(otus, correlation)) +
  geom_boxplot(fill = "#75aadb", alpha = 0.6, width = 0.4, outlier.shape = NA, position = position_nudge(x = -0.12)) +
  geom_jitter(aes(color = group), size = 0.5, width = 0.18, height = 0) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_color_manual(name = "", values = c("#a1b4c3", "#fbc47e")) +
  theme_classic() + theme(legend.position = "top")
ggplot2::ggsave(file.path(out_dir, "plot.png"), p041, height = 5, width = 10, dpi = 150, bg = "white")

# 已知回退：
# gghalves 的 transformation= 在 ggplot2 4 下失败，整理版用箱线+jitter 近似。

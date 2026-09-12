# =============================================================================
# 批量箱线显著性图
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
library(ggpubr)
library(ggsignif)
library(rstatix)

# 构造模拟数据：
data <- data.frame(
  TIME_IA = runif(10, min = 0.05, max = 0.4),
  TIME_ISM = runif(10, min = -0.2, max = 0.1),
  TIME_ISS = runif(10, min = 0.1, max = 0.5),
  TIME_IE = runif(10, min = -0.25, max = 0.2),
  TIME_IR = runif(10, min = 0.2, max = 0.6)
)

# 长宽数据转换：
data_long <- pivot_longer(data, cols = everything(),
                          names_to = "group", values_to = "Score")

data_long$group <- factor(data_long$group, levels = colnames(data))
# 计算显著性：
# 批量t检验：
stat.test <- data_long %>%
  wilcox_test(
    Score ~ group,
    p.adjust.method = "bonferroni"
  )

# 绘图：
colors <- c('#eb4b3a', "#48bad0", "#1a9781",
            "#355783", "#ef9a80")
p <- ggplot(data_long)+
  # 箱线图：
  geom_boxplot(aes(group, Score, color = group))+
  # 抖动散点：
  geom_jitter(aes(group, Score, color = group), width = 0.01)+
  # 颜色模式：
  scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                "#355783", "#ef9a80"))+
  xlab("")+
  # 主题：
  theme_classic()+
  theme(legend.position = "none",
        # x轴字体、颜色、角度调整：
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                                   color = colors))

# 根据显著性检验结果，添加显著性标记：
x_value <- rep(1:4, 4:1)
y_value <- rep(apply(data, 2, max)[1:4], 4:1) + 0.01
y_value <- y_value + c(0.03*1:4, 0.03*1:3, 0.03*1:2, 0.03)
color_value <- c(colors[2:5], colors[3:5], colors[4:5], colors[5])

for (i in 1:nrow(stat.test)) {
  if (stat.test$p.adj.signif[i] != "ns") {
    y_tmp <- y_value[i]
    p <- p+annotate(geom = "text",
                    label = stat.test$p.adj.signif[i],
                    x = x_value[i],
                    y = y_tmp,
                    color = color_value[i])
  }
}
p

ggsave(file.path(out_dir, "single_plot.pdf"), height = 4, width = 4)


# 循环绘制多图：
p_list <- list()
for (j in 1:6) {
  # 构造模拟数据：
  data <- data.frame(
    TIME_IA = runif(10, min = 0.05, max = 0.4),
    TIME_ISM = runif(10, min = -0.2, max = 0.1),
    TIME_ISS = runif(10, min = 0.1, max = 0.5),
    TIME_IE = runif(10, min = -0.25, max = 0.2),
    TIME_IR = runif(10, min = 0.2, max = 0.6)
  )

  # 长宽数据转换：
  data_long <- pivot_longer(data, cols = everything(),
                            names_to = "group", values_to = "Score")

  data_long$group <- factor(data_long$group, levels = colnames(data))
  # 计算显著性：
  # 批量t检验：
  stat.test <- data_long %>%
    wilcox_test(
      Score ~ group,
      p.adjust.method = "bonferroni"
    )

  # 绘图：
  colors <- c('#eb4b3a', "#48bad0", "#1a9781",
              "#355783", "#ef9a80")
  p <- ggplot(data_long)+
    # 箱线图：
    geom_boxplot(aes(group, Score, color = group))+
    # 抖动散点：
    geom_jitter(aes(group, Score, color = group), width = 0.01)+
    # 颜色模式：
    scale_color_manual(values = c('#eb4b3a', "#48bad0", "#1a9781",
                                  "#355783", "#ef9a80"))+
    xlab("")+
    # 主题：
    theme_classic()+
    theme(legend.position = "none",
          # x轴字体、颜色、角度调整：
          axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                                     color = colors))

  # 根据显著性检验结果，添加显著性标记：
  x_value <- rep(1:4, 4:1)
  y_value <- rep(apply(data, 2, max)[1:4], 4:1) + 0.01
  y_value <- y_value + c(0.03*1:4, 0.03*1:3, 0.03*1:2, 0.03)
  color_value <- c(colors[2:5], colors[3:5], colors[4:5], colors[5])

  for (i in 1:nrow(stat.test)) {
    if (stat.test$p.adj.signif[i] != "ns") {
      y_tmp <- y_value[i]
      p <- p+annotate(geom = "text",
                      label = stat.test$p.adj.signif[i],
                      x = x_value[i],
                      y = y_tmp,
                      color = color_value[i])
    }
  }
  p_list[[j]] <- p
}

library(cowplot)

plot_grid(plotlist = p_list, ncol = 3)

ggsave(file.path(out_dir, "all_plot.pdf"), height = 5.5, width = 9)

# =============================================================================
# 双向误差棒代谢标记散点
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
# 数据路径：file.path(root, "data", ...)
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

library(ggplot2)

############ 模拟数据构建 ----------------
# x值：
set.seed(123)
xmin <- runif(40)
xmax <- xmin + runif(40, 0.01, 0.08)

# y值：
set.seed(234)
ymin <- runif(40)
ymax <- ymin + runif(40, 0.01, 0.08)

# 数据框：
data <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
data <- data %>%
  mutate(xmean = (xmin+xmax)/2, ymean = (ymin+ymax)/2)

# 分组信息：
data$group1 <- rep(c(LETTERS[1:5], LETTERS[1:5]), 4)
data$group2 <- rep(c("E10", "E12", "E15", "E18"), each = 10)
data$group3 <- rep(c(rep("WT", 5), rep("AK", 5)), 4)

################# 绘图 --------------
ggplot(data, aes(x = xmean, y = ymean)) +
  # 横向误差棒：
  geom_errorbarh(aes(xmin = xmin, xmax = xmax,
                     color = group1, alpha = group2),
                 height = 0.01, linewidth = 0.2) +
  # 纵向误差棒：
  geom_errorbar(aes(ymin = ymin, ymax = ymax,
                    color = group1, alpha = group2),
                width = 0.01, linewidth = 0.2) +
  # 散点:
  geom_point(aes(fill = group1, alpha = group2, shape = group3))+
  # 斜线：
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  # 颜色模式：
  scale_fill_manual(values = c("#ffffff", "#ffae47", "#1e78c2", "#d44966", "#4b9258"))+
  scale_color_manual(values = c("#ffffff", "#ffae47", "#1e78c2", "#d44966", "#4b9258"))+
  scale_shape_manual(values = c(21, 22))+
  scale_alpha_manual(values = c(0.4, 0.6, 0.8, 1))+
  labs(x = "Fraction labeled aKG",
       y = "Fraction labeled glutamate") +
  theme_classic()+
  theme(legend.position = "none")

ragg::agg_png(file.path(out_dir, "xy_errorbar.png"), width = 4.2, height = 4, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "xy_errorbar.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "xy_errorbar.pdf"), last_plot(), width = 4.2, height = 4)











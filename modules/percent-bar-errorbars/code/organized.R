# =============================================================================
# 祖先成分百分比柱加误差棒
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」高分SCI图表复现 044百分比柱状图+误差棒
# 本地复现：drafts/ks-shengxin-huitu-repro/044-percent-bar
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
# 数据路径：file.path(root, "data", ...)
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############ 百分比柱状图+误差棒 #############

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

# 构建模拟数据：
data <- data.frame(
  sample = paste0("sample", 1:30),
  A = sample(40:100, 30, replace = T),
  B = sample(30:80, 30, replace = T),
  C = sample(10:20, 30, replace = T))

# 这里的数值假定为均值：
data_long <- data %>%
  pivot_longer(cols = !sample,
               names_to = "group",
               values_to = "mean_value")

# 随机构建标准差数据，当然你也可以认为是四分位间距：
data_long$sd <- sample(1:5, 90, replace = T)
data_long$max_value <- data_long$mean_value + data_long$sd
data_long$min_value <- data_long$mean_value - data_long$sd

head(data_long)
# # A tibble: 6 × 6
#   sample  group mean_value    sd max_value min_value
#   <chr>   <chr>      <int> <int>     <int>     <int>
# 1 sample1 A             88     5        93        83
# 2 sample1 B             62     5        67        57
# 3 sample1 C              5     3         8         2
# 4 sample2 A             69     2        71        67
# 5 sample2 B             47     3        50        44
# 6 sample2 C              6     2         8         4


########### 绘图 #############
data_long$group <- factor(data_long$group, levels = c(LETTERS[3:1]))


ggplot(data_long,aes(sample, mean_value, fill = group))+
  # 绘制百分比柱状图：
  geom_bar(aes(sample, mean_value, fill = group),
           stat = "identity", position = "fill", color = "grey")+
  # 修改颜色模式：
  scale_fill_manual(values = c("#979797", "#feda77", "#6bafd7"))+
  # 修改坐标轴标签：
  ylab("Propotion of ancestry")+
  xlab("")+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank())+
  # 翻转坐标轴：
  coord_flip()


########### 误差棒添加 -------------
# 这里比较有意思的就是误差棒不能直接用position = fill去修改，
# 你们可以自己尝试一下，我们需要先进行百分比转换；
for (i in 1:30) {
  for (j in 2:3) {
    data_long$max_value[(i-1)*3+j] <- sum(data_long$mean_value[((i-1)*3+1):((i-1)*3+j-1)])+data_long$max_value[(i-1)*3+j]
    data_long$min_value[(i-1)*3+j] <- sum(data_long$mean_value[((i-1)*3+1):((i-1)*3+j-1)])+data_long$min_value[(i-1)*3+j]
  }
}

data_long <- data_long %>%
  group_by(sample) %>%
  mutate(max_propotion = max_value/sum(mean_value),
         min_propotion = min_value/sum(mean_value))



ggplot(data_long,aes(sample, fill = group))+
  # 绘制百分比柱状图：
  geom_bar(aes(sample, mean_value, fill = group), width = 0.8,
           stat = "identity", position = "fill", color = "#303030")+
  geom_errorbar(aes(ymin = min_propotion, ymax = max_propotion),
                color = "#303030", width = 0.2)+
  # 修改颜色模式：
  scale_fill_manual(values = c("#979797", "#feda77", "#6bafd7"))+
  # 修改坐标轴标签：
  ylab("Propotion of ancestry")+
  xlab("")+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank())+
  # 翻转坐标轴：
  coord_flip()

ragg::agg_png(file.path(out_dir, "ancestry_bar.png"), width = 5, height = 10, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "ancestry_bar.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "ancestry_bar.pdf"), last_plot(), width = 5, height = 10)











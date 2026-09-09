# =============================================================================
# 癌症转移生存甘特条
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」高分SCI图表复现 009甘特图
# 本地复现：drafts/ks-shengxin-huitu-repro/009-gantt
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

# 甘特图
library(ggplot2)
library(dplyr)

#################### 示例数据 #################
data <- read.csv(file.path(root, "data", "test_data.csv"))

data$means <- apply(data[,c(2,3)], 1, mean)

head(data)
#   CancerType  low high means metastasis number
# 1   Prostate 0.80 0.95 0.875    Primary    217
# 2   Prostate 0.60 0.80 0.700  Bone Met.     77
# 3   Prostate 0.45 0.70 0.575  Lung Met.     19
# 4   Prostate 0.25 0.80 0.525 Brain Met.      7
# 5   Prostate 0.27 0.50 0.385 Liver Met.     38
# 6       Lung 0.55 0.72 0.635    Primary    952

table(data$metastasis)

data$color <- gsub("Bone Met.", "#947559", data$metastasis)
data$color <- gsub("Lung Met.", "#8678b0", data$color)
data$color <- gsub("Brain Met.", "#9f436b", data$color)
data$color <- gsub("Liver Met.", "#9abc6d", data$color)
data$color[which(data$color == "Primary")] <- c("#c97d80","#bba7cb",
                                                "#c1d9ec", "#f1bac8")

#################### 绘图 #################
# 对组内数据，根据means从大到小重排：
data2 <- data %>% 
  arrange(CancerType, means) %>% 
  ungroup %>%
  mutate(id=rep(c(1:5),4))

# 设置颜色变量：
cols <- data$color
names(cols) <- cols


# 最终绘图代码：
ggplot(data2)+
  # 这行代码构建空坐标系用：
  geom_point(aes(x=means, y=CancerType), color = "white")+
  # 背景阴影1：
  geom_rect(aes(xmin=0, xmax=Inf, ymin=1.5, ymax=2.5), fill = "#e6e6e6")+
  # 背景阴影2：
  geom_rect(aes(xmin=0, xmax=Inf, ymin=3.5, ymax=4.5), fill = "#e6e6e6")+
  
  # 虚线：
  geom_linerange(aes(xmin = 0, xmax = high, 
                     y = CancerType, group = id),
                 position = position_dodge(width = 0.5), 
                 linetype = "dashed")+
  # 方块：
  geom_tile(aes(x = means, y = CancerType, 
                height = 0.4, width = high-low,
                group = id, fill = color),
            color = "black",
            position = position_dodge(width = 0.5), 
            size = 0.3)+
  # 均值点：
  geom_point(aes(x=means, y=CancerType, group = id, fill=color), 
             shape=21, color = "black",
             position = position_dodge(width = 0.5), 
             size=3)+
  # 标签：
  geom_text(aes(label = paste0(metastasis, "(",number,")"), 
                x=high+0.08, y=CancerType, group = id), 
            position = position_dodge(width = 0.5), size=3)+
  # 主题：
  theme_classic() +
  # 去掉图形与坐标轴间隙并设置x轴刻度：
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     breaks = seq(0, 1, 0.25),
                     limits = c(0, 1.15))+
  scale_fill_manual(values = cols)+
  # 修改x轴和y轴标签：
  xlab("Area Under Kaplan-Meier Plot of Overall Survival\n(40 months follow-up)")+
  ylab("")+
  # 自定义主题：
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(angle=90, hjust = 0.5, vjust = 2,
                                   size = 10, color = "black"))

ragg::agg_png(file.path(out_dir, "gantt_chart.png"), width = 7, height = 7, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "gantt_chart.png"), file.path(root, "preview.png"), overwrite = TRUE)






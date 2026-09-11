# =============================================================================
# 分组提琴显著性图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」045分组提琴图+显著性检验
# 本地复现：drafts/ks-shengxin-huitu-repro/045-grouped-violin-signif
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

# 载入R包
library(ggplot2)
# library(tidyverse) 已在文件头拆包

# 构建数据：
data <- data.frame(Forwarding = sample(1000:7000, 300, replace = T),
                   Buffering = sample(0:3000, 300, replace = T),
                   Reinforcing = sample(0:2500, 300, replace = T))

data$group <- rep(c("Stage-specific", "Tissue-specific(E15.5)", "Tissue-specific(P42)"),
                  each = 100)

data_long <- data %>% pivot_longer(-group,
                                   names_to = "group2",
                                   values_to = "value")

# 计算每种分组的上下四分位数，用于绘制中心的森林图：
data_long <- data_long %>%
  group_by(group, group2) %>%
  mutate(High = quantile(value, 0.75),
         Med = median(value),
         Low = quantile(value, 0.25))

data_long$group2 <- factor(data_long$group2, levels = c("Forwarding", "Buffering",
                                                        "Reinforcing"))

# 计算显著性：
library(rstatix)
library(ggpubr)

stat_t_test <- data_long %>%
  group_by(group) %>%
  t_test(value~group2, paired = T) %>%
  ungroup()

stat_t_test <- stat_t_test %>% add_xy_position(x = "group")
stat_t_test$y.position <- stat_t_test$y.position+1500

##### 绘图 ------------
ggplot(data_long, aes(group, value, fill = group2))+
  # trim参数可以调节提琴图的尾部：
  geom_violin(trim = FALSE)+
  # 中心的森林图：
  geom_linerange(aes(x = group, ymin = Low, ymax = High, group = group2),
                 position = position_dodge(width = 0.9))+
  geom_point(aes(x = group, y = Med, group = group2),
             position = position_dodge(width = 0.9), size = 3)+
  # 显著性：
  stat_pvalue_manual(
    stat_t_test, label = "p.adj",
    bracket.size = 0.3, # 粗细
    tip.length = 0.01  # 两边竖线的长度
  )+
  scale_x_discrete(expand = c(0,0))+
  annotate("rect", xmin = c(0.55, 1.55, 2.55), xmax = c(1.45, 2.45, 3.45),
           ymin = rep(12000, 3), ymax = rep(13000, 3), fill = "#f2f2f2",
           color = "black")+
  annotate("text", x = c(1:3), y = rep(12500, 3),
           label = c(expression(bolditalic("Stage-specific")),
                     expression(bolditalic("Tissue-specific(E15.5)")),
                     expression(bolditalic("Tissue-specific(P42)"))),
           color = "black", size = 3)+
  # 颜色：
  scale_fill_manual(name = "", values = c("#dd8653", "#59a5d7", "#aa65a4"))+
  # 主题：
  theme_classic()+
  theme(legend.position = "bottom",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 7)

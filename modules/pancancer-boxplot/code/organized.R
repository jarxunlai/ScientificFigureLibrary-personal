# =============================================================================
# 泛癌箱线显著性图
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
set.seed(20260912)
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
library(ggpubr)
# library(tidyverse) 已在文件头拆包
library(latex2exp)
library(rstatix)

# 独立模拟数据：演示跨组箱线图和检验标记，不是实际癌症队列。
data <- expand.grid(project=paste0("Cohort ", LETTERS[1:12]), group=c("Normal","Tumor"), replicate=seq_len(40))
data$expr <- rnorm(nrow(data), 15 + as.integer(data$project)/5 + (data$group=="Tumor")*rep(seq(-1.5,1.5,length.out=12),80), 1.2)

# 绘图：
ggplot(data)+
  # 基础图形：
  geom_boxplot(aes(project, expr, fill = group))+
  # 设置颜色：
  scale_fill_manual(values = c("Tumor" = "#d6503a", "Normal" = "#5488ef"))+
  # 设置主题：
  theme_bw()+
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 按中位数由高到低排列：
data_new <- data %>%
  group_by(project) %>%
  mutate(median = median(expr), group_max = max(expr)) %>%
  arrange(desc(median))

# 调整因子顺序：
data_new$project <- factor(data_new$project, levels = unique(data_new$project))

data_new1 <- data_new %>%
  group_by(project) %>%
  mutate(group_number = length(unique(group)))

# 此数据用于添加显著性检验的label：
stat.test <- data_new1[which(data_new1$group_number == 2),] %>%
  group_by(project) %>%
  pairwise_t_test(
    expr ~ group, paired = FALSE,
    p.adjust.method = "fdr"
  ) %>%
  add_xy_position(x = "project")

stat.test

# 重新画图：
ggplot()+
  # 基础图形：
  geom_boxplot(data = data_new, aes(x = project, y =expr, fill = group),outlier.shape = 21, outlier.fill = "white")+
  # 设置颜色：
  scale_fill_manual(name = NULL, values = c("Tumor" = "#d6503a", "Normal" = "#5488ef"))+
  # 设置主题：
  theme_bw()+
  labs(y = TeX("Simulated expression ($log_{2}$ scale)", bold = T))+
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.99, 0.99),
        legend.justification = c(1,1))+
  # 显著性检验：
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif",
    bracket.size = 0.3, # 粗细
    tip.length = 0.01  # 两边竖线的长度
  )

ggsave(file.path(out_dir, "box_plot.pdf"), height = 6, width = 10)

ggplot2::ggsave(file.path(root,"preview.png"), plot=ggplot2::last_plot(), width=10, height=6, dpi=200, bg="white")

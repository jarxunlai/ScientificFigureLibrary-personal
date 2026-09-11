# =============================================================================
# 复杂百分比柱状图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」013复杂百分比柱状图
# 本地复现：drafts/ks-shengxin-huitu-repro/013-complex-percent-bar
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

# 百分比柱状图：
library(ggplot2)
# library(tidyverse) 已在文件头拆包

# 构建数据：
data <- data.frame(Q4 = c(32, 18, 14), Q1 = c(16, 28, 20), 
                   group = paste0("group",1:3)) %>% 
  pivot_longer(cols = !group, names_to = "X", values_to = "count")

# 绘图：
p1 <- ggplot(data)+
  geom_bar(aes(rev(X), count, fill = group), color = "#f3f4f4",
           position = "fill", stat = "identity", size = 1)+
  # 修改填充颜色：
  scale_fill_manual(values = c("#f6a34a", "#f6ddb4", "#bbbdc0"),
                    # 图例标签：
                    labels=rev(c("TC 0(<1%)","TC 1(<5%)","TC 2+(>=5%)")))+
  # 添加星号注释：
  annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  annotate("text", x = 1.5, y = 1.05, label=expression("*"~italic("P=0.003")), size = 4)+
  # 标题和副标题：
  ggtitle("Immune cell\nPD-L1", subtitle = "(by SP142 IHC)")+
  # 难点：学会使用expression函数：
  scale_x_discrete(labels = c(expression(atop(bold("Q4"), "(n=74)")),
                              expression(atop(bold("Q1"), "(n=74)"))))+
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        panel.background = element_rect(fill = "#f3f4f4")
        )+
  # 图例调整：
  # 图例顺序：
  guides(fill=guide_legend(reverse=TRUE))+
  # 图例标题：
  labs(fill="IC level")

p1

p2 <- ggplot(data)+
  geom_bar(aes(rev(X), count, fill = group), color = "#f3f4f4",
           position = "fill", stat = "identity", size = 1)+
  # 修改填充颜色：
  scale_fill_manual(values = c("#925fa7", "#c5b3d1", "#bbbdc0"),
                    # 图例标签：
                    labels=rev(c("TC 0(<1%)","TC 1(<5%)","TC 2+(>=5%)")))+
  # 添加星号注释：
  annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  annotate("text", x = 1.5, y = 1.05, label=expression("*"~italic("P=0.003")), size = 4)+
  # 标题和副标题：
  ggtitle("Immune cell\nPD-L1", subtitle = "(by SP142 IHC)")+
  # 难点：学会使用expression函数：
  scale_x_discrete(labels = c(expression(atop(bold("Q4"), "(n=74)")),
                              expression(atop(bold("Q1"), "(n=74)"))))+
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        panel.background = element_rect(fill = "#f3f4f4")
  )+
  # 图例调整：
  # 图例顺序：
  guides(fill=guide_legend(reverse=TRUE))+
  # 图例标题：
  labs(fill="TC level")


p3 <- ggplot(data)+
  geom_bar(aes(rev(X), count, fill = group), color = "#f3f4f4",
           position = "fill", stat = "identity", size = 1)+
  # 修改填充颜色：
  scale_fill_manual(values = c("#00acd5", "#7b8bc3", "#352e6d"),
                    # 图例标签：
                    labels=rev(c("Desert","Excluded","Inflamed")))+
  # 添加星号注释：
  annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  annotate("text", x = 1.5, y = 1.05, label=expression("*"~italic("P=0.003")), size = 4)+
  # 标题和副标题：
  ggtitle("Immune cell\nPD-L1", subtitle = "(by SP142 IHC)")+
  # 难点：学会使用expression函数：
  scale_x_discrete(labels = c(expression(atop(bold("Q4"), "(n=74)")),
                              expression(atop(bold("Q1"), "(n=74)"))))+
  # x轴和y轴标签
  xlab("")+
  ylab("")+
  # 设置主题：
  theme_bw()+
  theme(panel.grid = element_blank(),
        # 标题和副标题居中
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        # 修改背景色：
        panel.background = element_rect(fill = "#f3f4f4")
  )+
  # 图例调整：
  # 图例顺序：
  guides(fill=guide_legend(reverse=TRUE))+
  # 图例标题：
  labs(fill="Immune phenotype")


# 图形和图例的合并：
# 图例合并是难点！！！
library(cowplot)

p1a <- p1 + theme(legend.position = "none")
p2a <- p2 + theme(legend.position = "none")
p3a <- p3 + theme(legend.position = "none")

# 先合并无图例组：
p <- plot_grid(p1a, p2a, p3a, nrow = 1)

# 提取图例：
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)
legend3 <- get_legend(p3)

p_new <- plot_grid(p, plot_grid(legend1,legend2, legend3, ncol = 1, 
                                align = "v"),
                   rel_widths = c(3, 1))
# p_new

ggsave(file.path(out_dir, "plot.pdf"), plot = p_new, height = 5, width = 10)

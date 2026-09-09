# =============================================================================
# 微生物增温效应森林点距图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」高分SCI图表复现 033个性化森林图
# 本地复现：drafts/ks-shengxin-huitu-repro/033-forest
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

library(ggplot2)

# 构造数据：
min <- round(runif(28, -1, 0.7), 2)
sd <- round(runif(28, 0.2, 0.3), 2)
max <- min+sd
med <- (min+max)/2

data <- as.data.frame(cbind(min, med, max))
set.seed(33)
data$x <- factor(read.csv(file.path(root, "data", "rownames.csv"), header = FALSE)[,1],
                 levels = rev(read.csv(file.path(root, "data", "rownames.csv"), header = FALSE)[,1]))
data$group <- paste0("group", c(rep(1, 14), rep(2, 3),
                                rep(3, 3), rep(4, 5),
                                rep(5, 3)))
data$group_col <- c(rep("#e7a40e", 14), rep("#78bee5", 3),
                    rep("#1c6891", 3), rep("#a59d70", 5),
                    rep("#4f4a30", 3))

data$p <- c(rep("", 5), rep("*", 8),
            rep("**", 9), rep("***", 6))[sample(1:28)]

data$p_col <- ifelse(data$med > 0,
                     ifelse(data$p != "", "Postive effect(P<0.05)", "Postive effect(P>=0.05)"),
                     ifelse(data$p != "", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)"))

head(data)
#     min    med   max              x  group group_col   p                   p_col
# 1  0.17  0.320  0.47  Acidobacteria group1   #e7a40e   *  Postive effect(P<0.05)
# 2  0.33  0.430  0.53 Actinobacteria group1   #e7a40e     Postive effect(P>=0.05)
# 3 -0.75 -0.645 -0.54  Bacteroidetes group1   #e7a40e  **  Negtive effect(P<0.05)
# 4 -0.81 -0.695 -0.58     Chlamydiae group1   #e7a40e  **  Negtive effect(P<0.05)
# 5 -0.73 -0.625 -0.52    Chloroflexi group1   #e7a40e ***  Negtive effect(P<0.05)
# 6  0.03  0.150  0.27     Firmicutes group1   #e7a40e   *  Postive effect(P<0.05)


ggplot(data)+
  # 0轴竖线：
  geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条：
  geom_linerange(aes(x, ymin = min, ymax = max, color = p_col), show.legend = F)+
  # 散点：
  geom_point(aes(x, med, color = p_col)) +
  # 显著性：
  geom_text(aes(x, y = max + 0.05, label = p, color = p_col), show.legend = F)+
  # 颜色：
  scale_color_manual(name = "",
                     values = c("Postive effect(P<0.05)" = "#d55e00",
                                "Postive effect(P>=0.05)" = "#ffbd88",
                                "Negtive effect(P<0.05)" = "#0072b2",
                                "Negtive effect(P>=0.05)" = "#7acfff"))+
  # 背景色：
  annotate("rect",
           xmin = c(0.5,3.5,8.5,11.5,14.5),
           xmax = c(3.5,8.5,11.5,14.5,28.5),
           ymin = -1, ymax = 1, alpha = 0.2, fill = rev(unique(data$group_col))) +
  # 调整x轴拓宽：
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("Warming effect size")+
  theme_bw()+
  theme(axis.text.y = element_text(color = rev(data$group_col)))+
  coord_flip()

ragg::agg_png(file.path(out_dir, "forest_effect.png"), width = 6, height = 6, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "forest_effect.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "forest_effect.pdf"), last_plot(), width = 6, height = 6)


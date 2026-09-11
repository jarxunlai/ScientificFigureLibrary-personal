# =============================================================================
# 富集气泡图进阶
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」14富集分析气泡图进阶
# 本地复现：drafts/ks-shengxin-huitu-repro/c14-gobubble
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

############### 数据处理 ############## 
# 载入R包： 
library(GOplot)

# 示例数据：
data(EC)
class(EC)
dim(EC$david)
dim(EC$genelist)

# 将上述两个数据结合：
circ <- circle_dat(EC$david,EC$genelist)
class(circ)

# 查看数据形式：
head(circ)

# category         ID              term count  genes      logFC adj_pval
# 1       BP GO:0007507 heart development    54   DLC1 -0.9707875 2.17e-06
# 2       BP GO:0007507 heart development    54   NRP2 -1.5153173 2.17e-06
# 3       BP GO:0007507 heart development    54   NRP1 -1.1412315 2.17e-06
# 4       BP GO:0007507 heart development    54   EDN1  1.3813006 2.17e-06
# 5       BP GO:0007507 heart development    54 PDLIM3 -0.8876939 2.17e-06
# 6       BP GO:0007507 heart development    54   GJA1 -0.8179480 2.17e-06
# zscore
# 1 -0.8164966
# 2 -0.8164966
# 3 -0.8164966
# 4 -0.8164966
# 5 -0.8164966
# 6 -0.8164966

# labels:设置标签显示的阈值。阈值是指-log(adj_Pvalue) -- 默认值=5
# 简单的理解就是：看y轴，这个值以上的圈圈保留：
pdf(file.path(out_dir, "GOBubble1.pdf"), height = 8, width = 10)
GOBubble(circ, labels = 4)
dev.off()

############### 分面气泡图 ############## 
# 略调整参数之后可以对图的布局、颜色等进行调整：
pdf(file.path(out_dir, "GOBubble2.pdf"), height = 10, width = 15)
GOBubble(circ, title = 'Bubble plot',  # 标题
         colour = c('#e16d38', '#99c355', '#4cace9'),  # 气泡颜色
         display = 'multiple', # 表示分面
         labels = 3)
dev.off()

# 背景着色：
pdf(file.path(out_dir, "GOBubble3.pdf"), height = 10, width = 15)
GOBubble(circ, title = 'Bubble plot',  # 标题
         colour = c('#e16d38', '#99c355', '#4cace9'),  # 气泡颜色
         display = 'multiple', # 表示分面
         labels = 3,
         bg.col = T  # 背景着色，按照点的颜色分配背景色
         )
dev.off()

############## 进阶版：用ggplot2绘制类似的气泡图 ##################
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# 剔除一些重复的GO ID：
circ2<-circ[!duplicated(circ$ID),-5]

ggplot(circ2, aes(x = zscore, y = -log10(adj_pval))) +
  geom_point(aes(size = count,color=category),alpha = 0.6) +
  scale_size(range = c(1,12)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  geom_text_repel(
    data = circ2[-log10(circ2$adj_pval) > 3,],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE )
ggsave(file.path(out_dir, "GOBubble_ggplot1.pdf"), height = 7, width = 7)

############## 分面 ##################
ggplot(circ2, aes(x = zscore, y = -log10(adj_pval))) +
  geom_point(aes(size = count, color=category), alpha = 0.6) +
  scale_size(range =c (1, 12)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = c("none")) +
  geom_text_repel(
    data = circ2[-log10(circ2$adj_pval) > 3,],
    aes(label = ID),
    size = 3,
    segment.color = "black", show.legend = FALSE ) +
  facet_grid(.~category)

ggsave(file.path(out_dir, "GOBubble_ggplot2.pdf"), height = 5, width = 8)

# 已知回退：
# 环境无 GOplot，用模拟 GO 表画分面气泡。

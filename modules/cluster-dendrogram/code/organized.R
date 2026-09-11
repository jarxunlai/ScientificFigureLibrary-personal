# =============================================================================
# 聚类树
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」13聚类树
# 本地复现：drafts/ks-shengxin-huitu-repro/c13-cluster-tree
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

################## 基础树状图 ######################
# 载入R包
# library(tidyverse) 已在文件头拆包

# 数据
head(mtcars)

# 使用三个变量来聚类：
mtcars %>% 
  select(mpg, cyl, disp) %>% 
  dist() %>%  # 计算距离矩阵
  hclust() %>%   # 根据距离矩阵聚类
  as.dendrogram() -> dend  # 构建图对象

# Plot
pdf(file.path(out_dir, "clust1.pdf"), height = 7, width = 7)
par(mar=c(7,3,1,1))  # Increase bottom margin to have the complete label
plot(dend)
dev.off()


################## 修改树的特定部分属性 ######################
library(dendextend)

# 修改分支和标签属性：
pdf(file.path(out_dir, "clust2.pdf"), height = 7, width = 7)
dend %>% 
  # 自定义分支：
  set("branches_col", "grey") %>%  # 颜色；
  set("branches_lwd", 3) %>%  # 粗细;
  # 自定义标签：
  set("labels_col", "orange") %>%  # 颜色；
  set("labels_cex", 0.8) %>%  # 标签大小;
  plot()
dev.off()

# 修改节点属性：
pdf(file.path(out_dir, "clust3.pdf"), height = 7, width = 7)
dend %>% 
  set("nodes_pch", 19)  %>%  # 节点的形状
  set("nodes_cex", 0.7) %>%  # 节点的大小
  set("nodes_col", "orange") %>%  # 节点的颜色
  plot()
dev.off()

# 只修改子叶节点属性：
pdf(file.path(out_dir, "clust4.pdf"), height = 7, width = 7)
dend %>% 
  set("leaves_pch", 19)  %>% 
  set("leaves_cex", 0.7) %>% 
  set("leaves_col", "skyblue") %>% 
  plot()
dev.off()


################## 突出显示cluster ######################
pdf(file.path(out_dir, "clust5.pdf"), height = 7, width = 7)
par(mar=c(1,1,1,7))
dend %>%
  # k设置分几个cluster，k等于几，颜色就要填几个；
  set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz=TRUE, axes=FALSE)  # horiz表示水平绘制，默认为False
# 画一竖线：v=200表示x轴坐标，lty=2表示虚线；
abline(v = 200, lty = 2)
dev.off()

# 修改一个cluster的背景色 -- rect.dendrogram函数
pdf(file.path(out_dir, "clust6.pdf"), height = 7, width = 7)
par(mar=c(9,1,1,1))
dend %>%
  set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(axes=FALSE)
rect.dendrogram(dend, k=3, 
                lty = 5, # 描边线条类型
                lwd = 0.1, # 描边粗细
                x=1, 
                col=rgb(0.1, 0.2, 0.4, 0.1) # 填充颜色
                ) 
dev.off()


################## 比较两个聚类树 ######################
# 用两种不同的聚类方法做2个树状图
d1 <- USArrests %>% dist() %>% hclust( method="average" ) %>% as.dendrogram()
d2 <- USArrests %>% dist() %>% hclust( method="complete" ) %>% as.dendrogram()

# 自定义这些聚类树的属性，并将它们放在一个列表中
dl <- dendlist(
  d1 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  d2 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
)

# 放在一起绘制：
pdf(file.path(out_dir, "clust7.pdf"), height = 7, width = 7)
tanglegram(dl, 
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=7,
           lwd=2
)
dev.off()

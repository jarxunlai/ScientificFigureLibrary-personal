# =============================================================================
# 曼哈顿图
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

# install.packages("CMplot") # 安装包，如果已经安装，此行可忽略。
library(CMplot)

data(pig60K)   # calculated p-values by MLM
head(pig60K)

# SNP Chromosome Position    trait1     trait2     trait3
# 1 ALGA0000009          1    52297 0.7738187 0.51194318 0.51194318
# 2 ALGA0000014          1    79763 0.7738187 0.51194318 0.51194318
# 3 ALGA0000021          1   209568 0.7583016 0.98405289 0.98405289
# 4 ALGA0000022          1   292758 0.7200305 0.48887140 0.48887140
# 5 ALGA0000046          1   747831 0.9736840 0.22096836 0.22096836
# 6 ALGA0000047          1   761957 0.9174565 0.05753712 0.05753712


################### 普通曼哈顿图 ###################
# 此时会出三张图：分别对应trait1、2、3
CMplot(pig60K,  # 数据
       plot.type="m",  # type="m"，绘制曼哈顿图
       LOG10=TRUE,  # p值取-log10
       threshold=c(1e-6,1e-4),  # 设置阈值并添加阈值线
       # 设置阈值线的线型和粗细和颜色：
       threshold.lty=c(1,2), 
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"),
       amplify=TRUE, # 是否放大显著的点
       signal.col=c("red","green"),  # 设置显著点的颜色
       signal.cex=c(1,1),  # 显著点的大小
       signal.pch=c(19,19),  # 显著点的形状
       #chr.den.col=NULL,  # 设置图的颜色 
       file.output=TRUE,  # 是否输出图片
       file="jpg", # 输出图片格式（pdf、jpg或tiff）
       memo="",  # 输出图片的名称
       dpi=300,  # 分辨率
       verbose=TRUE, # 是否输出log信息
       ylim=NULL  # 可以设置y轴显示范围
       )  


# 在曼哈顿图底部加上SNP密度图：
CMplot(pig60K,  # 数据
       plot.type="m",  # type="m"，绘制曼哈顿图
       LOG10=TRUE,  # p值取-log10
       threshold=c(1e-6,1e-4),  # 设置阈值并添加阈值线
       # 设置阈值线的线型和粗细和颜色：
       threshold.lty=c(1,2), 
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"),
       amplify=TRUE, # 是否放大显著的点
       signal.col=c("red","green"),  # 设置显著点的颜色
       signal.cex=c(1,1),  # 显著点的大小
       signal.pch=c(19,19),  # 显著点的形状
       # 添加SNP密度图：下面两行关键代码
       bin.size=1e6, # SNP密度图的窗口大小
       chr.den.col=c("darkgreen", "white", "red"),  # 设置SNP密度图的颜色 
       
       file.output=TRUE,  # 是否输出图片
       file="jpg", # 输出图片格式（pdf、jpg或tiff）
       memo="SNP",  # 输出图片的名称
       dpi=300,  # 分辨率
       verbose=TRUE, # 是否输出log信息
       ylim=NULL  # 可以设置y轴显示范围
)  


################### Multi-track 曼哈顿图 ###################
# 会生成两张图：一张是三个trait共用一个坐标轴；另一张是将三张图拼接起来：
CMplot(pig60K, 
       plot.type="m", 
       multracks=TRUE,  # 关键参数，其它参数和上面的一样
       threshold=c(1e-6,1e-4),
       threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), 
       threshold.col=c("black","grey"), 
       amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "white", "red"), 
       signal.col=c("red","green"),
       signal.cex=c(1,1),
       file="jpg",
       memo="multi_track",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE)


################### 环形曼哈顿图 ###################
CMplot(pig60K,
       plot.type="c",  # 关键参数：c表示绘制环形曼哈顿图
       chr.labels=paste("Chr",c(1:18,"X","Y"),sep=""),  # 染色体名称
       r=0.4,  # 圆形半径
       cir.legend=TRUE,  # 图例
       outward=FALSE,  # 点的朝向是否向外
       cir.legend.col="black",  # 图例颜色
       cir.chr.h=1.3,  # 染色体边界的高度
       chr.den.col="black",
       file="jpg",
       memo="circle",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE)


# 添加显著性线，及修改显著点的特征：
CMplot(pig60K,
       plot.type="c",  # 关键参数：c表示绘制环形曼哈顿图
       chr.labels=paste("Chr",c(1:18,"X","Y"),sep=""),  # 染色体名称
       r=0.4,  # 圆形半径
       cir.legend=TRUE,  # 图例
       threshold=c(1e-6,1e-4),
       amplify=TRUE,
       # 设置显著的点特征：
       threshold.lty=c(1,2),
       threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),
       # 添加SNP密度图
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,
       outward=FALSE,  # 点的朝向是否向外
       cir.legend.col="black",  # 图例颜色
       cir.chr.h=1.5,  # 染色体边界的高度
       file="jpg",
       memo="circle_2",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE)

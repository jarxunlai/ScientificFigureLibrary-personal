# =============================================================================
# 哑铃图
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

library(ggplot2)


##################### 读取数据 #####################
data <- read.csv(file.path(root, "data", "school_earnings.csv"))

# 查看数据
head(data)
#       School Women Men Gap
# 1       MIT    94 152  58
# 2  Stanford    96 151  55
# 3   Harvard   112 165  53
# 4    U.Penn    92 141  49
# 5 Princeton    90 137  47
# 6   Chicago    78 118  40


##################### 绘图 #####################
# 基础绘图：
ggplot(data, aes(x=Women, xend=Men, y=School))+
  # 用 ggplot2 线段和端点构成哑铃：
  geom_segment(aes(xend=Men, yend=School), colour="gray", linewidth=.5)+
  geom_point(aes(x=Women), colour="#b1cb41", size=2)+
  geom_point(aes(x=Men), colour="#45a0e2", size=2)+
  # 主题：
  theme_light()+
  # 去掉x轴细网格线：
  theme(panel.grid = element_blank()) +
  xlab("Annual Salary (in thousands)")

ggsave(file.path(out_dir, "dunbbell_plot1.png"), height = 7, width = 5)

# 修饰和美化：
ggplot(data, aes(x=Women, xend=Men, y=School))+
  # 用 ggplot2 线段和端点构成哑铃：
  geom_segment(aes(xend=Men, yend=School), colour="gray", linewidth=.5)+
  geom_point(aes(x=Women), colour="#b1cb41", size=2)+
  geom_point(aes(x=Men), colour="#45a0e2", size=2)+
  # 加上两个散点做外环，并让其随着size大小而变化：
  geom_point(aes(Women, School,size=Women),
             alpha=0.5,color="#b1cb41")+
  geom_point(aes(Men, School,size=Men),
             alpha=0.5,color="#45a0e2")+
  # 主题：
  theme_light()+
  # 去掉x轴细网格线：
  theme(panel.grid = element_blank(),
        # 去掉图例：
        legend.position = "none") +
  xlab("Annual Salary (in thousands)")

ggsave(file.path(out_dir, "dunbbell_plot2.png"), height = 7, width = 5)


# 如果你熟练掌握了ggplot2各种修饰方法，那么你可以自定义出任何你想要的样式：
ggplot(data, aes(x=Women, xend=Men, y=School))+
  # 用 ggplot2 线段和端点构成哑铃：
  geom_segment(aes(xend=Men, yend=School), colour="gray", linewidth=.5)+
  geom_point(aes(x=Women), colour="#b1cb41", size=2)+
  geom_point(aes(x=Men), colour="#45a0e2", size=2)+
  # 加上两个散点做外环，并让其随着size大小而变化：
  geom_point(aes(Women, School,size=Women),
             alpha=0.5,color="#b1cb41")+
  geom_point(aes(Men, School,size=Men),
             alpha=0.5,color="#45a0e2")+
  # 主题：
  theme_light()+
  # 修改网格线：
  theme(panel.grid =element_blank(),
        # 去掉图例：
        legend.position = "none") +
  xlab("Annual Salary (in thousands)")

ggsave(file.path(out_dir, "dunbbell_plot3.png"), height = 7, width = 6)

# 已知回退：
# 原生 ggalt DLL 在当前 Pixi R 会崩溃，用 ggplot segment+point 回退。

ggplot2::ggsave(file.path(root, "preview.png"), plot = ggplot2::last_plot(), width = 6, height = 7, dpi = 200, bg = "white")

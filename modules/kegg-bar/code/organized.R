# =============================================================================
# KEGG柱状图
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

library(ggplot2)

############################ 加载数据 ##########################
data <- read.csv(file.path(root, "data", "data.csv"))

# 数据非常简单，就是两列，一列是通路，一列是通路中的基因数量：
# 通路的类别也包含在第一列
head(data)
#                           Pathway Number
# 1             Cellular processes     NA
# 2       Transport and catabolism    945
# 3 Cellular community-prokaryotes     17
# 4  Cellular community-eukaryotes    493
# 5                  Cell motility    221
# 6          Cell growth and death    414


############################# 绘图 #############################
# 首先绘制出基本的雏形，后面在调整细节：
# [local-repro skipped dangling ggplot+] ggplot(data)+
  # 柱状图：注意stat参数
  geom_bar(aes(Pathway, Number), stat = "identity")+
  # 反转x和y
  coord_flip()

# 这个图的难点就在于如何调整细节：
# 加颜色：对数据进行分组：
na_index <- which(is.na(data$Number))

group <- c(NA, rep(data$Pathway[na_index[1]], na_index[2] - na_index[1] - 1),
           NA, rep(data$Pathway[na_index[2]], na_index[3] - na_index[2] - 1),
           NA, rep(data$Pathway[na_index[3]], na_index[4] - na_index[3] - 1),
           NA, rep(data$Pathway[na_index[4]], na_index[5] - na_index[4] - 1),
           NA, rep(data$Pathway[na_index[5]], na_index[6] - na_index[5] - 1),
           NA, rep(data$Pathway[na_index[6]], nrow(data) - na_index[6]))

data$Group <- group
data$Pathway <- factor(data$Pathway, levels = rev(data$Pathway))

table(data$Group)

colors <- c("black", rep("#9dd1c9", na_index[2] - na_index[1] - 1),
            "black", rep("#f2b06f", na_index[3] - na_index[2] - 1),
            "black", rep("#bebbd7", na_index[4] - na_index[3] - 1),
            "black", rep("#eb8776", na_index[5] - na_index[4] - 1),
            "black", rep("#88afcf", na_index[6] - na_index[5] - 1),
            "black", rep("#f4b76e", nrow(data) - na_index[6]))

face <- c("bold", rep(NULL, na_index[2] - na_index[1] - 1),
          "bold", rep(NULL, na_index[3] - na_index[2] - 1),
          "bold", rep(NULL, na_index[4] - na_index[3] - 1),
          "bold", rep(NULL, na_index[5] - na_index[4] - 1),
          "bold", rep(NULL, na_index[6] - na_index[5] - 1),
          "bold", rep(NULL, nrow(data) - na_index[6]))

ggplot(data, aes(Pathway, Number))+
  # 柱状图：注意stat参数
  geom_bar(aes(fill = Group), stat = "identity")+
  # 加数据标签：
  geom_text(aes(label = Number, y = Number + 50), size = 2) +
  # 设置颜色：
  scale_fill_manual(values = c("#9dd1c9","#f2b06f","#bebbd7",
                               "#eb8776","#88afcf","#f4b76e"))+
  # 反转x和y
  coord_flip()+
  # 设置主题:
  theme_bw()+
  # 设置坐标轴刻度：
  scale_y_continuous(breaks=seq(0,1500, 500))+
  # 设置坐标轴标题：
  ylab("Number of Gene")+
  xlab("")+
  theme(legend.position = "none", # 去掉图例
        # 修改网格线：
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(linetype = "dashed"),
        panel.grid.major.y = element_line(linetype = "dashed"),
        # 去掉y轴刻度：
        axis.ticks.y = element_blank(),
        # y轴标签：
        axis.text.y = element_text(face = "bold",
                                   color = rev(colors),
                                   hjust = 0, # 左对齐
                                   size = 8,lineheight = 2),
        # 标题居中：
        plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(family = "Times")
        )+
  ggtitle("KEGG pathway annotation")

ggsave(file.path(out_dir, "KEGG.pdf"), height = 7, width = 7)
ggsave(file.path(out_dir, "KEGG.png"), height = 7, width = 7)

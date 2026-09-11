# =============================================================================
# 云雨图
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

# 加载依赖包，如未安装，请先安装，如上代码。
library(ggplot2)
library(gghalves)
# library(tidyverse) 已在文件头拆包

# Species因子化，后续作图需要
iris$Species <- iris$Species %>% factor(.,levels=unique(.))
head(iris)

# 绘制云雨图，需要用到gghalves的两个函数，geom_half_violin()和geom_half_dotplot()，分别绘制“云”和“雨
ggplot(iris , aes(x = Species, y = Sepal.Length, fill = Species))+
  # 绘制小提琴图的一半:
  geom_half_violin(aes(fill = Species),
                   position = position_nudge(x = .15, y = 0), # 偏移中心的距离；    
                   adjust=1.5, 
                   trim=FALSE, # 是否修缮提琴尾部
                   colour=NA, # 描边颜色
                   side = 'r') +
  # 绘制蜂窝图的一半：
  geom_half_dotplot(aes(fill = Species), 
                    method="histodot", stackdir="down", 
                    dotsize = 0.55, # 点的大小；
                    position=position_nudge(x = .1, y = 0), # 偏移中心的距离； 
                    binwidth = 0.15, # 相当于调整点的大小；
                    colour = NA)+
  coord_flip()

ggsave(file.path(out_dir, "cloudRainPlot1.pdf"), height = 7, width = 7)

################# 使用see包 #################
# 安装see包# 
# install.packages("see")
library(see)
library(ggthemes)
ggplot(iris, aes(x = Species,
                 y =Sepal.Length,
                 fill = Species)) + 
  geom_violindot(color = NA,  # 去除小提琴边框 
                 dots_color = NA, # 去除圆点边框
                 dots_size  = 0.5, # 圆点大小                
                 binwidth = 0.15, # 窗口bin的宽度                  
                 scale="area", # 云图以“面积”来缩放，还有"count" or "width"                    
                 trim = F # 是否过滤头尾数据，默认为TRUE
                 ) +                                        
  ggthemes::theme_base()+ # 换主题                                       
  coord_flip()+ # 旋转换坐标轴                                        
  scale_fill_material() # 修改填充颜色

ggsave(file.path(out_dir, "cloudRainPlot2.pdf"), height = 7, width = 7)

################ 来个分组的吧！顺便加个箱线图！ ###################3
data("iris")

# 模拟一列时间信息：
iris$time <- factor(sample(1:3, 150, replace = T), levels = c(1,2,3))

# 绘制云雨图，需要用到gghalves的两个函数，geom_half_violin()和geom_half_dotplot()，分别绘制“云”和“雨
ggplot(iris , aes(x = time, y = Sepal.Length, fill = Species))+
  # 绘制小提琴图的一半:
  geom_half_violin(aes(fill = Species),alpha=0.5,
                   position = position_nudge(x = .15, y = 0), # 偏移中心的距离；    
                   adjust=1.5, 
                   trim=FALSE, # 是否修缮提琴尾部
                   colour=NA, # 描边颜色
                   side = 'r') +
  geom_boxplot(aes(fill = Species),
               width=0.15,  # 调整箱线图宽度；
               alpha=0.5,
               position = position_dodge(0.2), # 调整箱线图之间的间距；
               outlier.shape = NA # 不显示离群值散点
               )+
  # 绘制蜂窝图的一半：
  geom_half_dotplot(aes(fill = Species), alpha=0.5,
                    method="histodot", stackdir="down", 
                    dotsize = 0.55, # 点的大小；
                    position=position_nudge(x = -0.15, y = 0), # 偏移中心的距离； 
                    binwidth = 0.15, # 相当于调整点的大小；
                    colour = NA)+
  coord_flip()+
  theme_bw()

ggsave(file.path(out_dir, "cloudRainPlot3.pdf"), height = 7, width = 7)

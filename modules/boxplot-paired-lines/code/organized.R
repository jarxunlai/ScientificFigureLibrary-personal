# =============================================================================
# 带连线箱线图
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

# 载入R包：
library(ggplot2)
library(RColorBrewer)
library(ggpubr) 

# 构建模拟数据：
set.seed(2000)
data <- data.frame(BAI2013 = rnorm(60),
                   class = rep(rep(letters[1:3], each=10),2),
                   treatment = rep(c("elevated","ambient"),each=30),
                   index=rep(seq(1,30),2))

head(data)
#      BAI2013 class treatment index
# 1  0.9001420     a  elevated     1
# 2 -1.1733458     a  elevated     2
# 3 -0.8974854     a  elevated     3
# 4 -1.4445014     a  elevated     4
# 5 -0.3310136     a  elevated     5
# 6 -2.9006290     a  elevated     6

# 设置颜色模式：
palette <- c(brewer.pal(7,"Set2")[c(4,5)])


###### 使用ggpaired函数 ###########
# ggpaired是ggpurb包中绘制配对箱线图的函数
ggpaired(data, x = "treatment", y = "BAI2013",
         fill = "treatment",  # fill指定分组变量：
         palette = palette,  # 颜色
         line.color = "grey50", # 散点连线颜色
         line.size = 0.15,  # 连线粗细
         point.size = 1.5,  # 散点大小
         width=0.6,  # 箱线图宽度
         facet.by = "class",  # 分面变量
         short.panel.labs = FALSE)+
  # 添加显著性检验：
  stat_compare_means(paired = TRUE)+
  # 设置主题：
  theme_minimal()+
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks=element_line(colour="black"))

ggsave(file.path(out_dir, "boxplot_line.pdf"), height = 5, width = 8)


############### 使用ggplot2绘制 ###########
library(reshape2)
library(ggforce)
library(dplyr)


# 构建贝塞尔曲线所需数据：
type <- as.character(unique(data$class))

df_bezier <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_bezier) <- c("index","treatment","class","value")

for (i in 1:length(type)){
  data0<-data[data$class==type[i],]
  
  data1<-split(data0,data0$treatment)
  
  data2<-data.frame(ambient=data1$ambient[,1],
                    elevated=data1$elevated[,1],
                    index=data1$ambient[,4])
  
  colnames(data2)<-c(1,2,"index")
  data2$'1.3'<-data2$'1'
  data2$'1.7'<-data2$'2'
  
  data3 <- melt(data2,id="index",variable.name ="treatment")
  data3$treatment <- as.numeric((as.character(data3$treatment)))
  
  data4<-arrange(data3,index,treatment) 
  data4$class<-type[i]
  
  df_bezier<-rbind(df_bezier,data4)
}


# 绘图：
# [local-repro skipped dangling ggplot+] ggplot(data)+
    # 箱线图：
    geom_boxplot(aes(x = factor(treatment), y = BAI2013,fill=factor(treatment)),
                 width=0.35,position = position_dodge(0),size=0.5,outlier.size = 0) + 
    # 散点：
    geom_point(aes(x = factor(treatment), y = BAI2013,fill=factor(treatment)),
               shape=21,colour="black",size=2)+
    # geom_bezier函数创建贝塞尔曲线：
    geom_bezier(data=df_bezier,
                aes(x= treatment, y = value, 
                    group = index,linetype = 'cubic'),
                size=0.25,colour="grey20") +
    # 设置填充颜色
    scale_fill_manual(values=brewer.pal(7,"Set2")[c(4,5)])+
    # 添加显著性检验：
    stat_compare_means(aes(x = factor(treatment), y = BAI2013),
                       paired = TRUE, method = "t.test")+
    # 分面：
    facet_grid(.~class, labeller = label_both)+
    # X轴y轴标题：
    labs(x="treatment",y="Value")+
    # 设置主题：
    theme_minimal()+
    theme(text = element_text(family = "serif"),  # 全部改成新罗马字体
        strip.background = element_rect(fill="#eaeae0", color = "#dcddcf"),  # 分面颜色
        strip.text = element_text(size=15,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid.minor=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks=element_line(colour="black"))

ggsave(file.path(out_dir, "boxplot_besier.pdf"), height = 5, width = 8)

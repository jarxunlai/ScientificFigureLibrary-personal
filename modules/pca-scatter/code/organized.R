# =============================================================================
# PCA图
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

############ prcomp函数 ################
# retx表示返回score，scale表示要标准化
iris.pca <- prcomp(iris[,-5],scale=T,rank=4,retx=T) #相关矩阵分解

summary(iris.pca) #方差解释度
iris.pca$sdev #特征值的开方
iris.pca$rotation #特征向量，回归系数
iris.pca$x #样本得分score


############ princomp函数 ################
library(stats)
# 默认方差矩阵(cor=F),改为cor=T则结果与prcomp相同
iris.pca<-princomp(iris[,-5],cor=T,scores=T) 
# 各主成份的SVD值以及相对方差
summary(iris.pca) 
# 特征向量，回归系数
iris.pca$loading 
iris.pca$score


############### 绘图 ####################
library(ggplot2)
library(ggrepel)

# 使用R自带iris数据集（150*5，行为样本，列为特征）
iris_input <- iris 

# 设置样本名
rownames(iris_input) <- paste("sample",1:nrow(iris_input),sep = "") 

# 查看数据集前几行
head(iris_input) 

# 进行PCA（scale. = TRUE表示分析前对数据进行归一化）
pca1 <- prcomp(iris_input[,-ncol(iris_input)],center = TRUE,scale. = TRUE)

df1 <- pca1$x # 提取PC score
df1 <- as.data.frame(df1) # 注意：如果不转成数据框形式后续绘图时会报错

summ1 <- summary(pca1)
summ1

# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

# 绘制PCA得分图
library(ggplot2)
ggplot(data = df1,aes(x = PC1,y = PC2,color = iris_input$Species))+
  # 绘制置信椭圆：
  stat_ellipse(aes(fill = iris_input$Species),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ 
  # 绘制散点：
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","orange","pink"))+
  #scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(file.path(out_dir, "pca01.pdf"),height = 5,width = 7)


################### 使用ggbiplot函数 #########################
library(ggbiplot)

# 注意：这个函数需要传入的是pca之后的结果list
ggbiplot(pca1, obs.scale = 1, var.scale = 1,
         # 分组信息：
         groups = iris_input$Species, 
         # 是否显示置信椭圆：
         ellipse = T, 
         # 是否显示中心的圆：
         circle = F) +
  scale_color_discrete(name = '') +

  theme_bw()+
  theme(legend.direction = 'horizontal', legend.position = 'top')

# 其实这个函数也是基于ggplot来做的，节约了代码！
# 而且ggplot的有些属性是直接可以加在上面的！
ggsave(file.path(out_dir, "pca02.pdf"),height = 5,width = 7)

############### 使用factoextra包 ######################
library(factoextra)

fviz_pca_biplot(res.pca, label = "var", 
                habillage=iris$Species,
                addEllipses=TRUE, 
                ellipse.level=0.95,
                ggtheme = theme_minimal())

ggsave(file.path(out_dir, "pca03.pdf"),height = 7,width = 10)

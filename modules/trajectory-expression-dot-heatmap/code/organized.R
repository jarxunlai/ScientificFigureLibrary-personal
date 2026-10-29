# 分轨迹绿色气泡表达热图 / Trajectory-faceted green expression dot heatmap
# 科学问题：不同轨迹阶段的基因在各细胞类型中的表达强度与阳性比例如何分布？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "trajectory-expression-dot-heatmap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复紧凑细胞类型网格、绿色气泡面积层级和三阶段分面间隙；输入/派生值均为合成示例，不冒充来源结果。
# 轨迹分面绿色气泡图：细胞类型 × 27 基因，面积编码阳性比例。
# 原作者脚本原样保留；本脚本只消费 data/expression_dotplot.csv，不运行单细胞上游。
library(ggplot2)
library(grid)


source(file.path(module_root,"code/export.R"))

## 1. 输入和原图顺序 ====
d<-read.csv(file.path(module_root,"data/expression_dotplot.csv"))
d$Gene<-factor(d$Gene,levels=rev(unique(d$Gene)))
d$CellType<-factor(d$CellType,levels=unique(d$CellType))
d$Trajectory<-factor(d$Trajectory,levels=c("Early","Intermediate","Late"))
write.csv(d,file.path(draft,"result/03_plot_input.csv"),row.names=FALSE)

## 2. 绘图：自由高度分面、外侧阶段标签、黑色细描边 ====
p<-ggplot(d,aes(CellType,Gene))+
  geom_point(aes(fill=Avg_Expr,size=Pct_Expr0),shape=21,colour="#333333",stroke=.18)+
  facet_grid(Trajectory~.,scales="free_y",space="free_y",switch="y")+
  scale_size_area(limits=c(0,1),breaks=c(.2,.4,.6,.8,1),max_size=4.2,labels=scales::percent)+
  scale_fill_distiller(palette="Greens",direction=1,limits=c(0,3.5),breaks=1:3)+
  scale_x_discrete(expand=expansion(add=.6))+
  labs(x=NULL,y=NULL,fill="Avg\nexpr",size="Pct nuclei\nexpr > 0")+
  theme_bw(base_size=8.5,base_family="Arial")+
  theme(panel.grid=element_blank(),panel.border=element_rect(colour="black",linewidth=.35),
        strip.placement="outside",strip.background.y=element_rect(fill="#d9d9d9",colour=NA),
        strip.text.y.left=element_text(angle=90,size=8,face="bold"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour="black",size=7.5),
        axis.text.y=element_text(size=8,colour="black"),axis.ticks=element_blank(),
        legend.title=element_text(size=8),legend.text=element_text(size=7.5),legend.key.height=unit(3,"mm"),
        panel.spacing.y=unit(1.5,"mm"),plot.margin=margin(2,2,2,2,"mm"))+
  guides(fill=guide_colourbar(barheight=unit(28,"mm"),barwidth=unit(3.5,"mm"),frame.colour="black",ticks=FALSE,order=1),
         size=guide_legend(order=2,override.aes=list(fill="white")))
save_gg(p,file.path(draft,"figure/03_trajectory_dotplot"),103,168)

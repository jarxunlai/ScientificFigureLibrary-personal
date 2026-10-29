# 三种表格森林图 / Three table-integrated forest plots
# 科学问题：效应量及其不确定性在疾病、人群或细胞群之间如何比较？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "table-integrated-forest")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复三类森林图表格字号/列宽及区间层级，保留截断箭头、双色HR和免疫亚群编号；输入/派生值均为合成示例，不冒充来源结果。
# 复现 3 种森林图｜纯 ggplot2 从零搭建
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1021)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/autoimmune.csv"))
df$y <- rev(seq_len(nrow(df)))
df$col <- ifelse(df$beta>.6,"#ed0e0e",ifelse(df$beta>.2,"#98c0e0","black"))
## 3. 森林图一：表格列 + 分区背景 + 截断区间箭头 ====
p <- ggplot(df,aes(y=y)) + annotate("rect",xmin=.2,xmax=.6,ymin=.5,ymax=21.5,fill="#e8f0f6") + annotate("rect",xmin=.6,xmax=1.6,ymin=.5,ymax=21.5,fill="#fdeaea") +
  geom_vline(xintercept=c(0,.2,.6),linetype=3,linewidth=.3,colour=c("grey30","#367ab1","red")) +
  geom_segment(aes(x=pmax(lower,-.8),xend=pmin(upper,1.6),yend=y),linewidth=.25,colour="grey45") +
  geom_segment(data=subset(df,lower < -0.8),aes(x=-.7,xend=-.8,yend=y),arrow=arrow(type="closed",length=unit(1,"mm")),linewidth=.25) +
  geom_segment(data=subset(df,upper>1.6),aes(x=1.5,xend=1.6,yend=y),arrow=arrow(type="closed",length=unit(1,"mm")),linewidth=.25) +
  geom_point(aes(x=beta,fill=col),shape=22,size=1.7,stroke=.3) + scale_fill_identity() +
  geom_text(aes(x=-5,label=Disease),hjust=0,size=2.5) + geom_text(aes(x=-3,label=Count),hjust=0,size=2.5) + geom_text(aes(x=-2.4,label=sprintf("%.2f (%.2f)",Mean,SD)),hjust=0,size=2.5) +
  geom_text(aes(x=2,label=sprintf("%.2f (%.2f, %.2f)",beta,lower,upper)),hjust=0,size=2.5) + geom_text(aes(x=3.9,label=format.pval(padj,digits=2,eps=.0001)),hjust=0,size=2.5) +
  annotate("text",x=c(-5,-3,-2.4,2,3.9),y=22.5,label=c("Disease","Count","Mean (SD)","beta (95% CI)","padj"),hjust=0,fontface="bold",size=2.8) +
  annotate("text",x=c(-.4,.4,1.1),y=22,label=c("Lower","Moderate","High"),size=2.6,colour=c("black","#367ab1","red")) +
  scale_x_continuous(limits=c(-5,4.8),breaks=seq(-.8,1.6,.4)) + scale_y_continuous(limits=c(.3,23),breaks=NULL) + labs(x="beta",y=NULL) +
  theme_void(base_family="Arial") + theme(axis.text.x=element_text(size=6),axis.title.x=element_text(size=7),axis.ticks.x=element_line(linewidth=.2),plot.margin=margin(4,4,4,4))
save_plot(p,file.path(fig_dir,"21a_autoimmune_forest_table"),215,135)
## 4. 森林图二：双人群HR与两列P值 ====
dd <- read.csv(file.path(module_root,"data/cancer_hr.csv"));dd$y <- match(dd$Cancer,rev(unique(dd$Cancer)))
dd$yp <- dd$y+ifelse(dd$Group=="DM(-); AI-IR(+)",.15,-.15)
stripes <- unique(dd[c("Cancer","y","Cases")])
p2 <- ggplot(dd,aes(y=yp)) + geom_rect(data=subset(stripes,y%%2==0),aes(ymin=y-.5,ymax=y+.5),xmin=-3,xmax=6.5,fill="grey94",inherit.aes=FALSE) +
  geom_vline(xintercept=1,linetype=3,linewidth=.3) + geom_segment(aes(x=lower,xend=upper,yend=yp),linewidth=.25,colour="grey50") +
  geom_point(aes(x=HR,fill=Group,colour=Group),shape=23,size=1.4,stroke=.35) + scale_fill_manual(values=c("DM(-); AI-IR(+)"="white","DM(+)"="black")) + scale_colour_manual(values=c("DM(-); AI-IR(+)"="#e46070","DM(+)"="black")) +
  geom_text(data=stripes,aes(x=-3,y=y,label=Cancer),inherit.aes=FALSE,hjust=0,size=2.5) + geom_text(data=stripes,aes(x=-.1,y=y,label=Cases),inherit.aes=FALSE,hjust=1,size=2.5) +
  geom_text(aes(x=ifelse(Group=="DM(-); AI-IR(+)",4.3,5.8),y=y,label=format.pval(Pvalue,digits=2,eps=.0001)),size=2.2) +
  annotate("text",x=c(-3,-.7,4.3,5.8),y=max(dd$y)+1,label=c("Cancer site","Cases","AI-IR (+) P","DM (+) P"),hjust=c(0,0,.5,.5),size=2.8,fontface="bold") +
  scale_x_continuous(limits=c(-3,6.5),breaks=seq(.5,3.5,.5)) + scale_y_continuous(limits=c(.3,max(dd$y)+1.5),breaks=NULL) + labs(x="Hazard ratio (95% CI)",y=NULL) +
  theme_void() + theme(axis.text.x=element_text(size=6),axis.title.x=element_text(size=7),legend.position="bottom",legend.text=element_text(size=6))
save_plot(p2,file.path(fig_dir,"21b_cancer_dual_hr"),210,155)
## 5. 森林图三：彩色编号免疫亚群与多变量协变量 ====
dd <- read.csv(file.path(module_root,"data/immune_hr.csv"));dd$y<-rev(seq_len(nrow(dd)))
pal <- c("Classical Monocyte"="#378c4e","Nonclassical Monocyte"="#7abb5e","FOLR2+ TAM"="#f4cccc","IFN-TAM"="#d8579a","Inflam-TAM"="#a599cb","RTM-TAM"="#7364aa","TREM2+ TAM"="#006cda","cDC1 (CLEC9A+)"="#b6daff","cDC2 (CD1C+)"="#e1a6cc","cDC3 (LAMP3+)"="#924800")
p3 <- ggplot(dd,aes(y=y)) + geom_rect(data=subset(dd,y%%2==0),aes(ymin=y-.5,ymax=y+.5),xmin=-1.5,xmax=4,fill="grey92",inherit.aes=FALSE) +
  geom_vline(xintercept=1,linetype=3,linewidth=.3) + geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.3) + geom_point(aes(x=HR),shape=15,size=1.5) +
  geom_text(aes(x=-.25,label=sprintf("%.2f (%.2f–%.2f)",HR,lower,upper)),size=2.3) +
  geom_point(data=subset(dd,Celltype!=""),aes(x=-1.35,fill=Celltype),shape=21,colour="transparent",size=3) + scale_fill_manual(values=pal) +
  geom_text(data=subset(dd,Celltype!=""),aes(x=-1.35,label=Number),colour="white",size=2.2,fontface="bold") +
  geom_text(data=subset(dd,Celltype==""),aes(x=-1.45,label=Variable),hjust=0,size=2.5) +
  geom_text(aes(x=2.6,label=format.pval(Pvalue,digits=2,eps=.0001)),hjust=0,size=2.5) +
  annotate("text",x=c(-1.45,-.25,2.6),y=max(dd$y)+1,label=c("Variable","HR (95% CI)","P value"),hjust=c(0,.5,0),fontface="bold",size=2.8) +
  scale_x_continuous(limits=c(-1.5,3.3),breaks=c(.5,1,1.5,2,2.5)) + scale_y_continuous(limits=c(.3,max(dd$y)+1.5),breaks=NULL) +
  labs(x="Hazard ratio",y=NULL) + theme_void() + theme(axis.text.x=element_text(size=6),legend.position="right",legend.text=element_text(size=6),legend.key.height=unit(3.5,"mm"))
save_plot(p3,file.path(fig_dir,"21c_immune_multivariable_forest"),210,108)

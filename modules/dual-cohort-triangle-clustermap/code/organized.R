# 双队列三角相关性聚类热图 / Dual-cohort triangular clustermap
# 科学问题：两个队列中细胞状态之间的相关结构是否一致？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "dual-cohort-triangle-clustermap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复双侧同序标签、上下树和双色标，修正整矩阵偏黄的相关结构；输入/派生值均为合成示例，不冒充来源结果。
Sys.setlocale("LC_CTYPE","English_United States.utf8")
sf<-sub("^--file=","",args[grep("^--file=",args)][1])

source(file.path(module_root,"code/export.R"))
library(circlize);library(RColorBrewer)
a<-read.csv(file.path(module_root,"data/cohort_GEO.csv"),row.names=1,check.names=FALSE)
b<-read.csv(file.path(module_root,"data/cohort_inhouse.csv"),row.names=1,check.names=FALSE)
# 原模拟表的共同偏移过强；正负载荷恢复跨模块的负相关区，不抄来源r值。
loading<-rep(c(1,1,-1,-1,1,-1),length.out=ncol(a))
a<-sweep(scale(a),2,loading,"*");b<-sweep(scale(b),2,loading,"*")
ra<-cor(a);rb<-cor(b);hc<-hclust(as.dist(1-(ra+rb)/2),method="average")
ord<-hc$order;ra<-ra[ord,ord];rb<-rb[ord,ord];n<-ncol(ra)
tree<-ggdendro::dendro_data(as.dendrogram(hc))$segments
upper<-colorRamp2(seq(-1,1,length.out=11),rev(brewer.pal(11,"Spectral")))
lower<-colorRamp2(seq(-1,1,length.out=11),rev(brewer.pal(11,"PuOr")))
write.csv(ra,file.path(entry,"result/04_visual_GEO_correlations.csv"))
write.csv(rb,file.path(entry,"result/04_visual_inhouse_correlations.csv"))
save_base({
 par(mar=c(0,0,0,0),family="Arial",ps=12,xaxs="i",yaxs="i")
 plot.new();plot.window(xlim=c(-13,n+14),ylim=c(-7,n+7),asp=1)
 for(i in 1:n)for(j in 1:n)rect(j-.5,n-i+.5,j+.5,n-i+1.5,
   col=if(i<j)upper(ra[i,j])else if(i>j)lower(rb[i,j])else "white",border=NA)
 for(sign in c(-1,1)){
   baseline<-if(sign==1)n+1.1 else -.1
   segments(tree$x,baseline+sign*tree$y/max(tree$y)*4,tree$xend,baseline+sign*tree$yend/max(tree$y)*4,lwd=.6,col="#444444")
 }
 text(.05,n:1,rownames(ra),adj=1,cex=.65);text(n+.95,n:1,rownames(ra),adj=0,cex=.65)
 text(.5,n+5.9,"In-house / GSE165897",adj=0,cex=.8,font=2)
 text(n,-5.8,"GSE165897",adj=1,cex=.8,font=2)
 z<-seq(-1,1,length.out=101)
 for(k in 1:2){y<-n+1+(k-1)*2.1;pal<-if(k==1)upper else lower
   rect(-12+head(z+1,-1)*3.2,y,-12+tail(z+1,-1)*3.2,y+.6,col=pal(head(z,-1)),border=NA)
   text(-12+(c(-1,0,1)+1)*3.2,y-.45,c(-1,0,1),cex=.5)
   text(-12,y+1,if(k==1)"GEO r"else"In-house r",adj=0,cex=.6)
 }
},file.path(entry,"figure/04_dual_cohort_triangle_clustermap"),205,180)

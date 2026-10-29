# 组织蛋白丰度内向环形图 / Inward circular tissue proteome abundance
# 科学问题：组织来源、蛋白丰度与定位评分如何在同一环形结构中展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "inward-circular-abundance-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复向内丰度柱、灰度GLS外轨、组织标签和图例之间的比例与字体层级；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 21组织蛋白图：分组宽度按蛋白数；灰度GLS；向内柱长=丰度。
## 采用circlize矢量矩形；不是点云，不需把柱形栅格化。
#############################################################################


source(file.path(module_root,"code/export.R"))
library(circlize)
d<-read.csv(file.path(module_root,"data/tissue_proteins.csv"));fig<-file.path(entry,"figure")
group_order<-names(sort(table(d$Group),decreasing=TRUE))
d$GLS_cat<-cut(d$GLS,breaks=c(0,2,3,4,Inf),labels=c("GLS1","GLS2","GLS3","GLS4"),right=FALSE,include.lowest=TRUE)
d<-d[order(match(d$Group,group_order),-as.integer(d$GLS_cat),-d$Abundance),]
colors<-setNames(c("#F65A3A","#7BC242","#4C63B6","#FFE31A","#C066B0","#73C7D9","#F28A1C","#7A5AAE","#95C93D","#E61E8C","#4DB4E6","#FFD100","#B7D92D","#F65A3A","#6CC36C","#D97A26","#6666CC","#F2A8C3","#28B7B2","#D95A4A","#B22222"),group_order)
gls_colors<-c(GLS1="#D9D9D9",GLS2="#999999",GLS3="#5A5A5A",GLS4="#040608")
counts<-as.numeric(table(factor(d$Group,levels=group_order)))
write.csv(d,file.path(entry,"result/11_tissue_proteins_sorted.csv"),row.names=FALSE)
write.csv(data.frame(Group=group_order,n_proteins=counts),file.path(entry,"result/11_tissue_counts.csv"),row.names=FALSE)

## 1. 三层轨道与组织图例，按组织内GLS再按丰度排列 ====
save_base({
  par(ps=11)
  circos.clear()
  circos.par(start.degree=90,gap.after=c(rep(.9,20),10),cell.padding=c(0,0,0,0),
             track.margin=c(.002,.002),points.overflow.warning=FALSE,
             canvas.xlim=c(-1.12,1.62),canvas.ylim=c(-1.13,1.16))
  circos.initialize(factors=factor(group_order,levels=group_order),
                    xlim=cbind(rep(0,21),counts))
  circos.track(ylim=c(0,1),track.height=.035,bg.col="white",bg.border="black",bg.lwd=.55,
    panel.fun=function(x,y){
      g<-CELL_META$sector.index;k<-match(g,group_order)
      circos.text(CELL_META$xcenter,if(k<=11)2.25 else 2.7,g,
                  facing=if(k<=11)"bending.inside"else"clockwise",niceFacing=TRUE,
                  adj=if(k<=11)c(.5,.5)else c(0,.5),cex=if(k<=11).76 else .60)
    })
  circos.track(ylim=c(0,1),track.height=.035,bg.border="black",bg.lwd=.55,
    panel.fun=function(x,y){
      z<-d[d$Group==CELL_META$sector.index,]
      for(i in seq_len(nrow(z)))circos.rect(i-1,0,i,1,col=gls_colors[as.character(z$GLS_cat[i])],border=NA)
    })
  circos.track(ylim=c(-1.04,0),track.height=.66,bg.border="black",bg.lwd=.5,
    panel.fun=function(x,y){
      g<-CELL_META$sector.index;z<-d[d$Group==g,]
      for(i in seq_len(nrow(z)))circos.rect(i-1+.015,-z$Abundance[i],i-.015,0,col=colors[g],border=NA)
      if(g==group_order[1]){
        circos.yaxis(side="left",at=seq(-1,0,.2),labels=as.character(seq(1,0,-.2)),labels.cex=5/7,tick.length=convert_x(1,"mm"))
      }
    })
  legend(1.19,.94,legend=group_order,fill=colors,border=NA,bty="n",cex=5/7,y.intersp=.82,x.intersp=.5)
  legend(1.19,-.68,legend=names(gls_colors)[4:1],fill=gls_colors[4:1],border=NA,bty="n",cex=5/7,y.intersp=.85,x.intersp=.5)
  circos.clear()
},file.path(fig,"11_inward_tissue_proteome"),180,170)

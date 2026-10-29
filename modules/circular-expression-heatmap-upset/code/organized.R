# 环形表达热图与中心 UpSet / Circular expression heatmap with an UpSet inset
# 科学问题：多个功能基因集合的交集组成与其表达模式如何联动展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-expression-heatmap-upset")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：扩大中心UpSet并恢复集合名称/计数层级，外圈仍用五组原生表达轨道；输入/派生值均为合成示例，不冒充来源结果。
# 分组环形表达热图 + 中心 UpSet。集合成员表为权威输入，交集大小由同一成员表计算。
# 热图同一基因允许跨组重复，表达值一致；UpSet 使用 unique 基因，不重复计数。
library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)


source(file.path(module_root,"code/export.R"))

## 1. 输入、成员矩阵与实际交集 ====
d<-read.csv(file.path(module_root,"data/grouped_expression.csv"))
membership<-read.csv(file.path(module_root,"data/set_membership.csv"),check.names=FALSE)
groups<-c("Cell_Cycle","Signaling","Hypoxia","Immune","Metabolism")
d$Group<-factor(d$Group,levels=groups)
d<-d[order(d$Group,as.integer(sub("Gene","",d$Gene))),]
m<-as.matrix(d[,paste0("Exp_",1:5)]);rownames(m)<-d$Gene
sets<-lapply(groups,function(g) membership$Gene[membership[,g]])
names(sets)<-groups
cm<-make_comb_mat(sets)
# 顺序按来源视觉从小交集到大交集递增，再按 degree，图中计数和保存表共用这一个顺序。
comb_order<-order(comb_size(cm),comb_degree(cm),decreasing=FALSE)
comb_table<-data.frame(Combination=comb_name(cm),Degree=comb_degree(cm),Count=comb_size(cm))
write.csv(comb_table[comb_order,],file.path(draft,"result/10_intersections.csv"),row.names=FALSE)
write.csv(d,file.path(draft,"result/10_ring_expression.csv"),row.names=FALSE)
col_fun<-colorRamp2(seq(-2,2,length.out=9),brewer.pal(9,"RdBu"))
gcol<-setNames(brewer.pal(9,"RdBu")[c(1,3,5,7,9)],groups)

## 2. 五组环形 + 固定中心视口内的标准 UpSet ====
save_draw({
  par(mar=c(0,0,0,0),family="Arial")
  circos.clear()
  circos.par(start.degree=60,gap.after=c(5,5,5,5,30),track.margin=c(0,.008),cell.padding=c(0,0,0,0),
             canvas.xlim=c(-1.24,1.38),canvas.ylim=c(-1.24,1.24),points.overflow.warning=FALSE)
  circos.heatmap(m,split=d$Group,cluster=FALSE,col=col_fun,track.height=.23,cell.border="white",cell.lwd=.4,
                 bg.border="black",bg.lwd=.5,rownames.side="outside",rownames.cex=.64)
  circos.track(track.index=get.current.track.index(),bg.border=NA,panel.fun=function(x,y){
    if(CELL_META$sector.numeric.index==5) {
      yy<-seq(CELL_META$cell.ylim[1],CELL_META$cell.ylim[2],length.out=6);yy<-(head(yy,-1)+tail(yy,-1))/2
      circos.text(rep(CELL_META$xlim[2],5)+convert_x(1.5,"mm"),yy,rev(colnames(m)),facing="inside",adj=c(0,.5),cex=.64)
    }
  })
  circos.track(ylim=c(0,1),track.height=.065,bg.col=adjustcolor(gcol,alpha.f=.32),bg.border="#555555",bg.lwd=.5,
                panel.fun=function(x,y) circos.text(CELL_META$xcenter,.5,CELL_META$sector.index,facing="bending.inside",niceFacing=TRUE,cex=.78))
  # 环中心图不是截图贴图，所有条形/连接点/标签均由成员矩阵原生渲染。
  pushViewport(viewport(x=.476,y=.52,width=.37,height=.26))
  upset<-UpSet(cm,set_order=groups,comb_order=comb_order,comb_col="#4393C3",pt_size=unit(.8,"mm"),lwd=.6,
               row_names_gp=gpar(fontsize=5.4,fontfamily="Helvetica"),
               top_annotation=upset_top_annotation(cm,gp=gpar(fill="#4393C3"),height=unit(12,"mm"),
                  annotation_name_gp=gpar(fontsize=5.3),axis_param=list(gp=gpar(fontsize=4))),
               left_annotation=upset_left_annotation(cm,gp=gpar(fill=adjustcolor(gcol,alpha.f=.35)),width=unit(12,"mm"),
                  annotation_name_gp=gpar(fontsize=5.3),axis_param=list(gp=gpar(fontsize=4))))
  draw(upset,newpage=FALSE,padding=unit(c(1,1,1,1),"mm"),background="transparent")
  popViewport()
  grid.text("Morphogenic\ngenes",x=.476,y=.36,gp=gpar(fontsize=9,fontface="bold",fontfamily="Arial"))
  lg<-Legend(title="Expression",col_fun=col_fun,at=c(-2,-1,0,1,2),legend_height=unit(25,"mm"),
             labels_gp=gpar(fontsize=6),title_gp=gpar(fontsize=6))
  draw(lg,x=unit(.965,"npc"),y=unit(.93,"npc"),just=c("right","top"))
  circos.clear()
},file.path(draft,"figure/10_circular_upset"),180,180)

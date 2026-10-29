# 上下三角双队列气泡相关性热图 / Dual-cohort triangular bubble correlation heatmaps
# 科学问题：M0 与 M1 两类样本的特征相关结构如何在同一坐标框架内比较？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "dual-cohort-bubble-triangle")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：为异色显著性版补回顶端和左侧三层注释，双队列两套图例保持独立；输入/派生值均为合成示例，不冒充来源结果。
# M0/M1 双三角气泡图：同色全矩阵 + 异色显著性筛选版。
# 原网页注释与循环方向不一致；本实现明确 i>j 为 M0、i<j 为 M1，相关基于真实输入观测。
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)


source(file.path(module_root,"code/export.R"))

## 1. Spearman 相关及逐对检验：P 值不是从 r 再算相关 ====
a<-read.csv(file.path(module_root,"data/M0_observations.csv"),check.names=FALSE)[,-1]
b<-read.csv(file.path(module_root,"data/M1_observations.csv"),check.names=FALSE)[,-1]
ann<-read.csv(file.path(module_root,"data/feature_annotations.csv"))
r0<-cor(a,method="spearman");r1<-cor(b,method="spearman")
n<-ncol(a);p0<-matrix(1,n,n,dimnames=dimnames(r0));p1<-p0
for(i in 1:n) for(j in i:n) {
  p0[i,j]<-p0[j,i]<-cor.test(a[,i],a[,j],method="spearman",exact=FALSE)$p.value
  p1[i,j]<-p1[j,i]<-cor.test(b[,i],b[,j],method="spearman",exact=FALSE)$p.value
}
m<-r0;m[upper.tri(m)]<-r1[upper.tri(r1)];diag(m)<-1
pm<-p0;pm[upper.tri(pm)]<-p1[upper.tri(p1)]
for(pair in list(c("M0_r", "r0"),c("M1_r","r1"),c("M0_p","p0"),c("M1_p","p1"))) write.csv(get(pair[2]),file.path(draft,"result",paste0("09_",pair[1],".csv")))
cat_col<-c(SR="#00bad2",Immune="#f06dac",Metabolism="#67c575",NE="#d8b152",Proliferation="#8674cc",Other="#c0b6a1")
tumor_col<-c(Prostate="#f7d56b","Pan-cancer"="#d18a56",Breast="#90a8b5",Lung="#8b6ca8",Melanoma="#ef9c9a",Other="#b5b5b5")
cycle_col<-c(Yes="#b96aaa",No="#a9d8bc")
same_col<-colorRamp2(seq(-1,1,length.out=11),brewer.pal(11,"RdBu"))
upper_col<-colorRamp2(c(-1,0,1),c("#75c500","white","#e43a21"))
lower_col<-colorRamp2(c(-1,0,1),c("#f6e213","white","#801b9c"))

## 2. 同色版：顶部类别条、左侧三层注释、完整相关强度 ====
ha<-HeatmapAnnotation(Category=ann$Category,col=list(Category=cat_col),show_annotation_name=FALSE,show_legend=FALSE,
                      gp=gpar(col="#cccccc",lwd=.3),simple_anno_size=unit(2.3,"mm"))
ra<-rowAnnotation(Category=ann$Category,TumorType=ann$TumorType,KEGGcycle=ann$KEGGcycle,
                  col=list(Category=cat_col,TumorType=tumor_col,KEGGcycle=cycle_col),
                  show_annotation_name=FALSE,simple_anno_size=unit(2.3,"mm"),gp=gpar(col="#eeeeee",lwd=.25),
                  annotation_legend_param=list(Category=list(nrow=2,title_gp=gpar(fontsize=6.8),labels_gp=gpar(fontsize=5)),
                    TumorType=list(nrow=2,title_gp=gpar(fontsize=6.8),labels_gp=gpar(fontsize=5)),KEGGcycle=list(nrow=2,title_gp=gpar(fontsize=6.8),labels_gp=gpar(fontsize=5))))
ht<-Heatmap(m,name="Spearman r",col=same_col,rect_gp=gpar(type="none"),
  top_annotation=ha,left_annotation=ra,cluster_rows=FALSE,cluster_columns=FALSE,
  row_names_side="left",column_names_side="top",row_names_gp=gpar(fontsize=6.4,fontfamily="Helvetica"),column_names_gp=gpar(fontsize=6.4,fontfamily="Helvetica"),
  cell_fun=function(j,i,x,y,width,height,fill){
    grid.rect(x,y,width,height,gp=gpar(fill=NA,col="#cccccc",lwd=.25))
    grid.circle(x,y,r=sqrt(abs(m[i,j]))*.43*min(unit.c(width,height)),gp=gpar(fill=same_col(m[i,j]),col=NA))
  },heatmap_legend_param=list(direction="horizontal",title_gp=gpar(fontsize=7.5),labels_gp=gpar(fontsize=6.8)))
save_draw({draw(ht,heatmap_legend_side="bottom",annotation_legend_side="bottom",merge_legends=TRUE,padding=unit(c(3,3,3,3),"mm"))},file.path(draft,"figure/09a_same_palette_annotated"),180,192)

## 3. 异色版：P<0.05 才显示气泡，缺失/不显著不映射为零 ====
ra2<-rowAnnotation(Category=ann$Category,TumorType=ann$TumorType,KEGGcycle=ann$KEGGcycle,
 col=list(Category=cat_col,TumorType=tumor_col,KEGGcycle=cycle_col),show_annotation_name=FALSE,show_legend=FALSE,
 simple_anno_size=unit(2.3,"mm"),gp=gpar(col="#eeeeee",lwd=.25))
ht2<-Heatmap(m,top_annotation=ha,left_annotation=ra2,rect_gp=gpar(type="none"),show_heatmap_legend=FALSE,cluster_rows=FALSE,cluster_columns=FALSE,
  row_names_side="left",column_names_side="top",row_names_gp=gpar(fontsize=6.4,fontfamily="Helvetica"),column_names_gp=gpar(fontsize=6.4,fontfamily="Helvetica"),
  cell_fun=function(j,i,x,y,width,height,fill){
    grid.rect(x,y,width,height,gp=gpar(fill=NA,col="#cccccc",lwd=.25))
    if(i==j) grid.circle(x,y,r=.43*min(unit.c(width,height)),gp=gpar(fill="#c4c4c4",col=NA))
    else if(pm[i,j]<.05) {
      col<-if(i>j) lower_col(m[i,j]) else upper_col(m[i,j])
      grid.circle(x,y,r=sqrt(abs(m[i,j]))*.43*min(unit.c(width,height)),gp=gpar(fill=col,col=NA))
    }
  })
lg1<-Legend(title="M1 (upper)",col_fun=upper_col,direction="horizontal",title_gp=gpar(fontsize=7.5),labels_gp=gpar(fontsize=6.8),legend_width=unit(38,"mm"))
lg0<-Legend(title="M0 (lower)",col_fun=lower_col,direction="horizontal",title_gp=gpar(fontsize=7.5),labels_gp=gpar(fontsize=6.8),legend_width=unit(38,"mm"))
save_draw({draw(ht2,annotation_legend_list=list(lg1,lg0),annotation_legend_side="bottom",padding=unit(c(3,3,3,3),"mm"))},file.path(draft,"figure/09b_dual_palette_significant"),180,187)

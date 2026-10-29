# 形态发生基因环状热图与集合重叠 / Circular morphogenic gene heatmap with set overlap
# 科学问题：三类动力学基因在五个心脏区域的表达与集合重叠如何对应？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-gene-heatmap-venn")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复环形热图、q菱形和组名层级的占比；保留两种来源配色和实际交集数；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 五区域环状表达热图、q值环与三集合Venn。
## 集合数从实际输入计算；中心Venn不作面积比例拟合，未安装eulerr。
#############################################################################


source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/gene_expression.csv"));fig<-file.path(entry,"figure")
groups<-c("Divergence","Curl","Acceleration")
regions<-c("Left_atrium","Right_atrium","Outflow_tract","Right_ventricle","Left_ventricle")
d<-d[order(match(d$Group,groups),d$qvalue),]
sets<-lapply(groups,function(g)unique(d$Gene[d$Group==g]));names(sets)<-groups
union_genes<-unique(d$Gene)
membership<-sapply(sets,function(s)union_genes%in%s)
pattern<-apply(membership,1,paste0,collapse="")
keys<-c("TRUEFALSEFALSE","FALSETRUEFALSE","FALSEFALSETRUE","TRUETRUEFALSE","TRUEFALSETRUE","FALSETRUETRUE","TRUETRUETRUE")
counts<-as.integer(table(factor(pattern,levels=keys)))
overlap<-data.frame(region=c("Divergence only","Curl only","Acceleration only","Divergence & Curl","Divergence & Acceleration","Curl & Acceleration","All three"),count=counts)
write.csv(overlap,file.path(entry,"result/04_gene_set_overlap.csv"),row.names=FALSE)
write.csv(d,file.path(entry,"result/04_expression_qvalues.csv"),row.names=FALSE)
sizes<-as.numeric(table(factor(d$Group,levels=groups)))
available<-360-3*6-24;spans<-sizes/sum(sizes)*available
qfun<-circlize::colorRamp2(c(0,.002,.005,.01),c("#FFF1E7","#F99C70","#97044D","#000000"))

## 1. 教程PiYG版与原图RdBu版分别保存，数据与布局完全共用 ====
for(palette in c("tutorial_piyg","source_rdbu")){
  expression_col<-if(palette=="tutorial_piyg")circlize::colorRamp2(c(0,.25,.5),c("#4DAC26","#F7F7F7","#D01C8B")) else circlize::colorRamp2(c(0,.25,.5),c("#2166AC","#F7F7F7","#B2182B"))
  group_col<-if(palette=="tutorial_piyg")c("#C3E3B2","#FFFFFF","#E8B6D3") else c("#8FB6DC","#F7F7F7","#F3B3A0")
  save_base({
  par(ps=11)
    canvas(c(-1.18,1.42),c(-1.18,1.18));a0<-60
    qxy<-matrix(numeric(0),ncol=2);qcolors<-character()
    for(k in seq_along(groups)){
      z<-d[d$Group==groups[k],];step<-spans[k]/nrow(z);a1<-a0-spans[k]
      for(i in seq_len(nrow(z))){
        a<-a0-(i-.5)*step
        for(j in seq_along(regions))sector(a-step/2,a+step/2,.72+(j-1)*.05,.72+j*.05,expression_col(z[[regions[j]]][i]),"white",.25)
        qxy<-rbind(qxy,xy(a,.685));qcolors<-c(qcolors,qfun(z$qvalue[i]))
        radial_text(a,.982,z$Gene[i],cex=5/7)
      }
      sector(a1,a0,.575,.635,group_col[k],"grey25",.5)
      arc_text((a0+a1)/2,.605,groups[k],cex=6/7)
      a0<-a1-6
    }
    raster_points(qxy,pch=18,col=qcolors,cex=.52)
    # 外圈的五行名置于顶部开口；图例标注与输入量纲一致。
    for(j in seq_along(regions)){
      p<-xy(90,.72+(j-.5)*.05);segments(p[1,1],p[1,2],p[1,1]+.03,p[1,2],lwd=.5)
      text(p[1,1]+.04,p[1,2],gsub("_"," ",regions[j]),adj=0,cex=5/7)
    }
    text(.04,.685,"q value",adj=0,cex=5/7)
    color_key(expression_col,c(0,.5),1.24,.67,height=.34,ticks=seq(0,.5,.1),title="Expression fraction")
    color_key(qfun,c(0,.01),1.24,-.95,height=.34,ticks=c(0,.005,.01),title="q value")
    # 当前拟态成员表的三个集合无交集：并列圆忠实显示实际集合，不伪造重叠。
    order_sets<-c(2,3,1);cc<-c(-.31,0,.31);rad<-c(.14,.11,.16)
    for(k in 1:3){a<-seq(0,2*pi,length.out=200)
      polygon(cc[k]+rad[k]*cos(a),.05+rad[k]*sin(a),col=adjustcolor(group_col[order_sets[k]],.65),border="grey20",lty=2,lwd=.65)
      text(cc[k],.05,length(sets[[order_sets[k]]]),cex=6/7)
    }
    text(0,-.405,"Morphogenic\ngenes",cex=1.1,font=2)
  },file.path(fig,paste0("04_morphogenic_",palette)),180,180)
}

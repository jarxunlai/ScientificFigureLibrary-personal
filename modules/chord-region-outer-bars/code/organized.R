# 作物区域和弦与管理措施外圈条形图 / Crop-region chord diagram with management bars
# 科学问题：作物、区域与驱动因子的关系如何与管理措施指标并置？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "chord-region-outer-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复内圈标签和带宽层级，补出独立内圈和弦版本；完整图仍保留外圈管理条形；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 内圈作物-地区与驱动因子和弦；外圈是独立的管理措施表。
## 使用同一几何坐标绘制，不把不匹配的类别强制当成同一扇区。
#############################################################################


source(file.path(module_root,"code/export.R"))
links<-read.csv(file.path(module_root,"data/inner_links.csv"));cr<-read.csv(file.path(module_root,"data/crop_regions.csv"))
outer<-read.csv(file.path(module_root,"data/outer_management.csv"));fig<-file.path(entry,"figure")
region_col<-setNames(c("#6BCCEF","#00D187","#009462","#FDC848","#9C9C9C","#951255"),unique(cr$Region))
crop_col<-c(Maize="#810084",Rice="#0075B6",Wheat="#000000")
nodes<-c(cr$node,unique(links$to));node_col<-setNames(c(region_col[cr$Region],rep("#525252",9)),nodes)
weight<-setNames(rep(0,length(nodes)),nodes)
for(i in seq_len(nrow(links))){weight[links$from[i]]<-weight[links$from[i]]+links$value[i];weight[links$to[i]]<-weight[links$to[i]]+links$value[i]}
span<-(360-length(nodes)*1.3)*weight/sum(weight)
start<-setNames(90-c(0,head(cumsum(span+1.3),-1)),nodes);end<-start-span
cursor<-start
links$a1<-links$a2<-links$da1<-links$da2<-NA_real_
for(i in seq_len(nrow(links))){
  f<-links$from[i];t<-links$to[i]
  df<-span[f]*links$value[i]/weight[f];dt<-span[t]*links$value[i]/weight[t]
  links$a1[i]<-cursor[f]-df/2;links$a2[i]<-cursor[t]-dt/2
  links$da1[i]<-df;links$da2[i]<-dt;cursor[f]<-cursor[f]-df;cursor[t]<-cursor[t]-dt
}
write.csv(links,file.path(entry,"result/06_chord_links.csv"),row.names=FALSE)
write.csv(outer,file.path(entry,"result/06_management_bars.csv"),row.names=FALSE)

## 1. 和弦带精确占据对应节点权重宽度 ====
for(view in c("full","inner")) save_base({
  par(ps=11)
  canvas(if(view=="full")c(-1.62,1.28)else c(-.87,.87),if(view=="full")c(-1.28,1.36)else c(-.87,.87))
  for(i in seq_len(nrow(links))){
    a1<-links$a1[i];a2<-links$a2[i];w1<-links$da1[i];w2<-links$da2[i]
    tt<-seq(0,1,length.out=45)
    p1<-xy(a1-w1/2,.66);p2<-xy(a2+w2/2,.66)
    p3<-xy(a2-w2/2,.66);p4<-xy(a1+w1/2,.66)
    b1<-cbind((1-tt)^3*p1[1]+tt^3*p2[1],(1-tt)^3*p1[2]+tt^3*p2[2])
    b2<-cbind((1-tt)^3*p3[1]+tt^3*p4[1],(1-tt)^3*p3[2]+tt^3*p4[2])
    aa<-seq(a2+w2/2,a2-w2/2,length.out=6)
    bb<-seq(a1+w1/2,a1-w1/2,length.out=6)
    polygon(rbind(b1,xy(aa,.66),b2,xy(bb,.66)),col=adjustcolor(node_col[links$from[i]],.34),border=NA)
  }
  for(k in seq_along(nodes)){
    sector(end[k],start[k],.67,.715,node_col[k],"white",.5)
    lab<-if(k<=nrow(cr))cr$Crop[k] else nodes[k]
    radial_text((start[k]+end[k])/2,.727,lab,cex=5/7)
  }
  # 外圈与和弦图不共享节点，依照Species独立分组和排序。
  if(view=="full") {
  step<-360/nrow(outer)
  for(i in seq_len(nrow(outer))){
    a<-90-(i-.5)*step;rr<-.87+outer$EF_mean[i]/max(outer$EF_mean)*.15
    sector(a-.45*step,a+.45*step,.87,rr,crop_col[outer$Species[i]],"white",.3)
    radial_text(a,rr+.012,outer$Type[i],cex=5/7)
    if(step>2.1)radial_text(a,(rr+.87)/2,sprintf("%.1f",outer$EF_mean[i]),"white",5/7)
  }
  legend(-1.56,1.26,legend=names(crop_col),fill=crop_col,bty="n",cex=5/7,title="Outer circle",title.adj=0)
  legend(-1.56,1.00,legend=c(names(region_col),"Driver"),fill=c(region_col,"#525252"),bty="n",cex=5/7,title="Inner circle",title.adj=0,y.intersp=.85)
  }
},file.path(fig,if(view=="full")"06_crop_region_chord_management"else"06_crop_region_chord"),180,180)

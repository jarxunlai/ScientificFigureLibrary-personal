# 抗氧化指标环形柱与误差线 / Circular antioxidant bars with error intervals
# 科学问题：四种处理下的氧化应激和抗氧化指标如何变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "polar-antioxidant-error-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：纠正四处理颜色语义：蓝、粉、浅蓝、橙；恢复分量纲扇区及误差线和字母；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 十一项抗氧化指标、四处理：独立量纲、误差线和显著性字母。
## 原文第七篇未提供代码，仅指向第二篇；本脚本独立实现全部图层。
#############################################################################


source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/antioxidant_traits.csv"));fig<-file.path(entry,"figure")
traits<-unique(d$category);groups<-c("CK","Cd","GSH","GSH+Cd")
bg<-c("#FFF9D2","#FBEE9D","#F9D546","#EBCBE0","#E49CBD","#D83D92","#E4EFE3","#C6DFC3","#4C9B31","#E2D5ED","#AA82CC")
colors<-setNames(c("#73C7E3","#E588D2","#C5E5F1","#F0AF41"),groups)
axis_max<-setNames(c(6,.4,50,800,400,800,120,20,20,12,5),traits)
span<-360/length(traits);gap<-5
write.csv(d,file.path(entry,"result/07_antioxidant_intervals.csv"),row.names=FALSE)

## 1. 误差来自明确se列；未执行上游显著性检验 ====
save_base({
  par(ps=11)
  canvas(c(-1.13,1.23),c(-1.13,1.18))
  for(k in seq_along(traits)) {
    z<-d[d$category==traits[k],];z<-z[match(groups,z$group),]
    a0<-90-(k-1)*span-gap/2;a1<-90-k*span+gap/2
    mx<-axis_max[traits[k]]
    sector(a1,a0,.24,.93,adjustcolor(bg[k],.33),"grey45",.45)
    sector(a1,a0,.98,1.016,bg[k],"grey35",.5)
    arc_text((a0+a1)/2,.954,paste0(traits[k]," (",z$unit[1],")"),cex=5/7)
    for(t in (1:3)/4){arc(a1,a0,.24+t*.69,"grey60",.35,3);p<-xy(a0-.4,.24+t*.69);text(p[1,1],p[1,2],format(mx*t,trim=TRUE),srt=a0-90,cex=5/7,adj=c(1,.5))}
    angles<-seq(a0-3.4,a1+3.4,length.out=4)
    for(j in 1:4){
      rr<-.24+z$value[j]/mx*.69;e<-z$se[j]/mx*.69
      sector(angles[j]-2.1,angles[j]+2.1,.24,rr,colors[z$group[j]],"grey30",.45)
      radial_segment(angles[j],rr-e,rr+e,lwd=.55)
      arc(angles[j]-1.1,angles[j]+1.1,rr+e,"black",.55)
      arc(angles[j]-1.1,angles[j]+1.1,rr-e,"black",.55)
      radial_text(angles[j],rr+e+.022,z$significance[j],cex=5/7,facing="tangent")
    }
  }
  legend(.90,1.12,legend=groups,fill=colors,border="grey40",bty="n",cex=5/7,y.intersp=.8)
},file.path(fig,"07_antioxidant_polar_errorbars"),180,180)

# 根系多指标极坐标条形图 / Polar bars for multiple root traits
# 科学问题：不同处理的根系指标在各自量纲内如何比较？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "polar-trait-error-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：补齐根系量纲标签，放大扇形主体、独立刻度与三处理显著性字母；输入/派生值均为合成示例，不冒充来源结果。
Sys.setlocale("LC_CTYPE","English_United States.utf8")
#############################################################################
## 七项根系指标：同扇区内共用量纲，跨扇区独立范围。
## 输入 data/root_traits.csv；字母为示例输入注释，不重算多重检验。
#############################################################################
script<-sub("^--file=","",args[grep("^--file=",args)][1])

source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/root_traits.csv")); fig<-file.path(entry,"figure")
traits<-c("Forks","Crossings","Length","SurfArea","AvgDiam","RootVolume","Tips")
bg<-c("#FFF6B8","#F9E5EC","#F6ADCC","#DE3D8F","#7ABF43","#B7D9B0","#DBE7D5")
colors<-c(KB="#76C6E6",BC="#BBDDEB",BAc="#E4F1F5")
axis_max<-c(Forks=10000,Crossings=6000,Length=10000,SurfArea=1500,AvgDiam=1,RootVolume=180,Tips=1800)
span<-360/length(traits); gap<-6
write.csv(d,file.path(entry,"result/02_root_traits.csv"),row.names=FALSE)

## 1. 分别输出原教程基本版和误差线版 ====
for(with_error in c(FALSE,TRUE)) {
  save_base({
  par(ps=11)
    canvas(c(-1.13,1.17),c(-1.13,1.17))
    for(k in seq_along(traits)) {
      z<-d[d$category==traits[k],];z<-z[match(names(colors),z$group),]
      a0<-90-(k-1)*span-gap/2; a1<-90-k*span+gap/2
      max_val<-axis_max[traits[k]]
      sector(a1,a0,.26,.94,adjustcolor(bg[k],alpha.f=.33),"grey50",.4)
      sector(a1,a0,.98,1.012,bg[k],"grey30",.5)
      label<-c(Forks="Forks",Crossings="Crossings",Length="Length(cm)",SurfArea="SurfArea(cm²)",AvgDiam="AvgDiam(mm)",RootVolume="RootVolume(cm³)",Tips="Tips")[traits[k]]
      arc_text((a0+a1)/2,.952,label,cex=.95)
      for(t in (1:3)/4){arc(a1,a0,.26+t*.68,"grey60",.35,3);p<-xy(a0-.7,.26+t*.68);text(p[1,1],p[1,2],format(max_val*t,trim=TRUE),srt=a0-90,cex=5/7,adj=c(1,.5))}
      angles<-seq(a0-6,a1+6,length.out=3)
      for(j in 1:3) {
        rr<-.26+z$value[j]/max_val*.68
        sector(angles[j]-4.2,angles[j]+4.2,.26,rr,colors[z$group[j]],"grey30",.5)
        if(with_error) {
          e<-z$se[j]/max_val*.68;radial_segment(angles[j],rr-e,rr+e,lwd=.6)
          arc(angles[j]-1.6,angles[j]+1.6,rr+e,"black",.6)
          arc(angles[j]-1.6,angles[j]+1.6,rr-e,"black",.6)
        }
        lab_r<-rr+(if(with_error)z$se[j]/max_val*.68 else 0)+.028
        radial_text(angles[j],lab_r,z$significance[j],cex=6/7,facing="tangent")
      }
    }
    legend(.85,1.12,legend=names(colors),fill=colors,border="grey40",bty="n",cex=5/7,xpd=NA,y.intersp=.85)
  },file.path(fig,if(with_error)"02_root_traits_errorbars" else "02_root_traits_polar"),180,180)
}

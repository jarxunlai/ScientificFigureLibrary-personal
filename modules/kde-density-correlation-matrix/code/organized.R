# 核密度相关矩阵R版本 / KDE correlation matrix in R
# 科学问题：多变量关联的方向与平滑二维密度如何一起查看？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "kde-density-correlation-matrix")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：修正七变量符号、统一0–5.5密度坐标和完整竖直相关色标，收紧单元间隙；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 七变量数据；Spearman数值从同一表计算 ====
d <- read.csv(file.path(dat,"nucleosome_features.csv"),check.names=FALSE)
vars <- names(d); n <- ncol(d); r <- cor(d,method="spearman")
labels <- expression(sp^{4*"+"},spd^{3*"+"},CoH^{3*"+"},PEG,Ca^{2*"+"},HP1*alpha,HP1*beta*"+SUV")
## 2. 两套调色保持原文；下三角为二维高斯核密度 ====
for(style in c("rainbow","red")) {
  tile_cols <- if(style=="rainbow") colorRampPalette(c("#00007f","blue","cyan","green","yellow","#ff7f00","red","#7f0000"))(256) else colorRampPalette(c("#fff5f0","#fcae91","#fb6a4a","#cb181d"))(256)
  dens_cols <- colorRampPalette(c("white","#17becf","#1f77b4","#2ca02c","yellow","#ff7f0e","#d62728","#8b0000"))(256)
  save_base({
    layout(cbind(matrix(1:(n*n),n,n,byrow=TRUE),rep(n*n+1,n)),widths=c(rep(1,n),.65))
    par(mar=c(.85,.85,.30,.20),mgp=c(.8,.18,0),tcl=-.15,ps=12,cex=.92,family="Arial")
    for(i in 1:n) for(j in 1:n) {
      if(i==j) {
        plot.new();text(.5,.5,labels[i],cex=1.1)
      } else if(i<j) {
        plot.new();val <- r[i,j]; idx <- round((pmax(.1,pmin(.8,val))-.1)/.7*255)+1
        rect(0,0,1,1,col=tile_cols[idx],border="grey40");text(.5,.5,sprintf("%.2f",val),cex=1.0,font=2)
      } else {
        xx <- d[[j]]; yy <- d[[i]]
        den <- MASS::kde2d(xx,yy,n=65)
        image(den$x,den$y,den$z,col=dens_cols,xlim=c(0,5.5),ylim=c(0,5.5),axes=FALSE,xlab="",ylab="",useRaster=TRUE,zlim=c(0,max(den$z)))
        box(col="grey45",lwd=.6)
        if(i==n) axis(1,at=c(0,2.5,5),cex.axis=.8)
        if(j==1) axis(2,at=c(0,2.5,5),cex.axis=.8)
      }
    }
    par(mar=c(5,0,5,1.4),ps=12,cex=.9);plot.new();plot.window(xlim=c(0,1),ylim=c(.1,.8),xaxs="i",yaxs="i")
    yy<-seq(.1,.8,length.out=257);rect(.10,head(yy,-1),.62,tail(yy,-1),col=tile_cols,border=NA)
    axis(4,at=seq(.1,.8,.1),cex.axis=.7);mtext("Spearman correlation",4,line=1.2,cex=.8)
  },file.path(fig,paste0("12_density_matrix_",style)),190,180)
}

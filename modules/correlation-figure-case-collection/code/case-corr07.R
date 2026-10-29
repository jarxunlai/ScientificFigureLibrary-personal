# 十种相关性图复刻案例总览 / Correlation case-study collection
# 科学问题：不同变量维度与数据结构应选用何种相关性表达？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "correlation-figure-case-collection")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：重建矩阵和相关网络的统一坐标，网络不再跨过下三角散点面板；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure/case-corr07")
dat <- file.path(module_root,"data/case-corr07")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 叶片性状、森林类型与叶习性 ====
d <- read.csv(file.path(dat,"leaf_traits.csv")); vars <- c("LeafN","Chl","SLA","StomatalDensity","PoreLength")
names_pretty <- c("Leaf N\nconcentration","Chl\nconcentration","SLA","Stomatal\ndensity","Pore length")
cols <- c(CB="#2d7977",DL="#95dddc",HH="#f4a98c",XSB="#cd0000")
shapes <- c(Deciduous=16,Evergreen=17)
## 2. 下三角散点矩阵，右上空位放置相关网络 ====
save_base({
 par(mar=c(0,0,0,0),family="Arial",ps=11)
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 # 网络严格位于矩阵上三角；目标节点接在对角单元的右上角，连线不穿过单元。
 par(fig=c(.07,.82,.10,.88),new=TRUE,mar=c(0,0,0,0),xaxs="i",yaxs="i")
 plot.new();plot.window(xlim=c(.5,5.5),ylim=c(.5,5.5))
 sx<-c(2.0,3.8,5.35);sy<-c(5.35,4.65,3.1)
 for(i in 1:3)for(j in 1:5){
  ct<-cor.test(d[[c("Gwmax","Vcmax","WUEi")[i]]],d[[vars[j]]]);t<-seq(0,1,length.out=80)
  xx<-(1-t)*sx[i]+t*(j+.48);yy<-(1-t)*sy[i]+t*(6-j+.48)+.14*sin(pi*t)
  lines(xx,yy,col=if(ct$estimate>0)"#2e4fdb"else"#fe0108",lty=if(ct$p.value<.05)1 else 3,lwd=.5+1.8*abs(ct$estimate))
 }
 points((1:5)+.48,(5:1)+.48,pch=21,bg="grey80",cex=.65)
 points(sx,sy,pch=21,bg="grey85",cex=.8);text(sx,sy-.17,c("Gwmax","Vcmax","WUEi"),cex=.85,font=2)
 for(i in 1:5)for(j in 1:i){
  par(fig=c(.07+(j-1)*.15,.07+j*.15,.10+(5-i)*.156,.10+(6-i)*.156),new=TRUE,
      mar=c(.75,.75,.15,.10),mgp=c(.7,.2,0),tcl=-.13,ps=11,cex=.7)
  if(i==j){plot.new();box(col="grey50");text(.5,.5,names_pretty[i],cex=.95)}else{
   plot(d[[vars[j]]],d[[vars[i]]],col=cols[d$Site],pch=shapes[d$LeafType],cex=.65,xlab="",ylab="",axes=FALSE)
   abline(lm(d[[vars[i]]]~d[[vars[j]]]),lwd=.7)
   if(i==5)axis(1,cex.axis=.75);if(j==1)axis(2,cex.axis=.75);box(col="grey50")
  }
 }
 par(fig=c(.84,1,.14,.88),new=TRUE,mar=c(0,0,0,0),ps=10,cex=1)
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
 legend(0,1,names(cols),col=cols,pch=16,bty="n",cex=.7)
 legend(0,.78,names(shapes),pch=shapes,bty="n",cex=.65)
 legend(0,.55,c("Pearson r > 0","Pearson r < 0","P < 0.05","P >= 0.05"),
   col=c("#2e4fdb","#fe0108","grey30","grey30"),lty=c(1,1,1,3),bty="n",cex=.65)
},file.path(fig,"07_grouped_scatter_matrix_network"),215,180)

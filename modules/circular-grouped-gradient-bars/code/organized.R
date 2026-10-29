# 分组筛选与连续渐变环形条形图 / Grouped screening and continuous circular bar plots
# 科学问题：分组筛选结果与连续变量排序分别适合怎样的环形编码？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-grouped-gradient-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复分组筛选的刻度开口和内外向柱形，连续城市图保持同一零环及数值配色；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 两个独立场景：分组筛选抑制率、连续城市指标。
## 所有负值由0参考环向内绘制；城市值及名称来自明确示例CSV。
#############################################################################


source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/screening.csv"));cities<-read.csv(file.path(module_root,"data/city_indicator.csv"))
fig<-file.path(entry,"figure")
groups<-unique(d$Group);colors<-setNames(c("#4798AC","#B5333B"),groups)
d$angle<-NA_real_;step<-170/36
d<-d[order(match(d$Group,groups),ifelse(d$Group==groups[1],-d$Value,d$Value)),]
for(k in 1:2){ix<-which(d$Group==groups[k]);d$angle[ix]<-c(-90,90)[k]-(seq_along(ix)-.5)*step}
write.csv(d,file.path(entry,"result/09_grouped_screening.csv"),row.names=FALSE)

## 1. 分组抑制率；零刻度与数值范围不因视觉效果改变 ====
save_base({
  par(ps=11)
  canvas(c(-1.22,1.22),c(-1.22,1.22))
  rzero<-.55;scale<-.006
  for(k in 1:2){
    ix<-which(d$Group==groups[k]);a0<-max(d$angle[ix])+step/2;a1<-min(d$angle[ix])-step/2
    arc(a1,a0,rzero,"grey30",.6)
    arc_text(mean(c(a0,a1)),.48,groups[k],6/7,colors[k])
    for(i in ix){
      rr<-rzero+d$Value[i]*scale
      sector(d$angle[i]-.43*step,d$angle[i]+.43*step,rzero,rr,colors[k],"white",.2)
      radial_text(d$angle[i],max(rzero,rr)+.02,d$Type[i],cex=5/7)
    }
  }
  radial_segment(90,rzero-20*scale,rzero+80*scale,lwd=.6)
  for(v in seq(-20,80,20)){y<-rzero+v*scale;segments(-.012,y,.012,y,lwd=.55);text(-.025,y,v,adj=1,cex=5/7)}
  text(.06,.98,"Inhibition (%)",srt=90,cex=6/7)
},file.path(fig,"09_grouped_screening_inhibition"),180,180)

## 2. 连续指标：所有柱从同一个0参考圆起，图例和柱共用color_fun ====
cities<-cities[order(cities$CPUG),]
cities$angle<-0-(seq_len(nrow(cities))-.5)*354/nrow(cities)
fun<-circlize::colorRamp2(seq(-2.3,2.3,length.out=11),RColorBrewer::brewer.pal(11,"RdYlBu"))
write.csv(cities,file.path(entry,"result/09_city_continuous_indicator.csv"),row.names=FALSE)
save_base({
  par(ps=11)
  canvas(c(-1.20,1.20),c(-1.20,1.20))
  rzero<-.67;scale<-.12;step<-354/nrow(cities)
  arc(0,-354,rzero,"grey30",.5)
  for(i in seq_len(nrow(cities))){
    rr<-rzero+cities$CPUG[i]*scale
    sector(cities$angle[i]-.46*step,cities$angle[i]+.46*step,rzero,rr,fun(cities$CPUG[i]),"white",.25)
    radial_text(cities$angle[i],max(rzero,rr)+.022,cities$City[i],cex=5/7)
  }
  color_key(fun,c(-2.3,2.3),-.04,-.34,width=.045,height=.68,ticks=-2:2,title="CPUG (K)")
  radial_segment(0,rzero-2*scale,rzero+2*scale,lwd=.55)
  for(v in -2:2){x<-rzero+v*scale;segments(x,-.012,x,.012,lwd=.5);text(x,.03,v,srt=90,cex=5/7)}
},file.path(fig,"09_continuous_city_bars"),180,180)

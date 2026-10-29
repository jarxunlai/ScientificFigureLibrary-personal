# 人群特征多轨环形图 / Cohort characteristics in a multitrack circular plot
# 科学问题：个体年龄、BMI和分类特征的分布如何在同一视图中对应？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "cohort-circular-multitrack")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：收紧主环画布，恢复270度环与开口摘要的视觉占比及样本标签；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 队列基线五轨环形图；输入 synthetic subjects.csv，输出独立 PNG/PDF。
## 年龄排序沿角度对应同一参与者；四分之一开口内展示真实输入的分布摘要。
#############################################################################



source(file.path(module_root,"code/export.R"))
d <- read.csv(file.path(module_root,"data/subjects.csv"))
d <- d[order(d$age),]
fig <- file.path(entry,"figure"); result <- file.path(entry,"result")
age_col <- "#631879"; bmi_col <- "#FF410D"
sex_col <- c(Female="#F2300F",Male="#E1BD6D")
eth_col <- c(Caucasian="#ECCBAE",Asian="#0B775E",Hispanics="#D69C4E",Black="#81A88D",Unknown="grey70")
iris_col <- c(IS="#F21A00",IR="#74A089",Unknown="grey70")
step <- 270/nrow(d); angle <- 90-(seq_len(nrow(d))-.5)*step
write.csv(d,file.path(result,"01_cohort_multitrack.csv"),row.names=FALSE)

## 1. 五层同心环与开口摘要；主图为一个协调的多轨视图 ====
save_base({
  par(ps=11)
  canvas(c(-1.15,1.17),c(-1.14,1.18))
  for(i in seq_len(nrow(d))) {
    a<-angle[i]; lo<-a-step/2; hi<-a+step/2
    sector(lo,hi,.81,1,"white","grey25",.35)
    radial_segment(a,.81,.81+(d$age[i]-20)/60*.19,age_col,.8)
    sector(lo,hi,.59,.79,"white","grey25",.35)
    radial_segment(a,.59,.59+(d$BMI[i]-16)/28*.20,bmi_col,.8)
    sector(lo,hi,.46,.57,sex_col[d$Gender[i]],"grey25",.35)
    sector(lo,hi,.33,.44,eth_col[d$Ethnicity[i]],"grey25",.35)
    sector(lo,hi,.20,.31,iris_col[d$IRIS[i]],"grey25",.35)
    radial_text(a,1.012,d$subject[i],cex=5/7)
  }
  # >100个点只栅格化点层；圆环、线、文字保留矢量。
  raster_points(xy(angle,.81+(d$age-20)/60*.19),col=age_col,cex=.36)
  raster_points(xy(angle,.59+(d$BMI-16)/28*.20),col=bmi_col,cex=.36)
  for(k in c(20,50,80))text(-.025,.81+(k-20)/60*.19,k,adj=1,cex=5/7)
  for(k in c(16,30,44))text(-.025,.59+(k-16)/28*.20,k,adj=1,cex=5/7)
  # 摘要栏与对应径向轨道对齐；箱线与堆积计数均直接来自d。
  for(v in c("age","BMI")) {
    x<-if(v=="age") -.88 else -.66
    col<-if(v=="age") age_col else bmi_col
    lim<-if(v=="age") c(20,80) else c(16,44)
    yy<-.03+(d[[v]]-lim[1])/diff(lim)*.66
    rect(x-.075,.03,x+.075,.69,border="grey70",lwd=.4)
    bin<-round(yy/.013); offsets<-ave(seq_along(bin),bin,FUN=function(z)(seq_along(z)-(length(z)+1)/2)*.012)
    raster_points(x+offsets,yy,pch=16,cex=.35,col=col)
    qs<-quantile(yy,c(0,.25,.5,.75,1)); segments(x,qs[1],x,qs[5],lwd=.5)
    rect(x-.06,qs[2],x+.06,qs[4],border="grey35",lwd=.5)
    segments(x-.06,qs[3],x+.06,qs[3],lwd=.9)
    text(x,.725,if(v=="age") "Age" else "BMI",cex=6/7)
  }
  for(k in 1:3) {
    x<-c(-.51,-.38,-.25)[k]
    variable<-c("Gender","Ethnicity","IRIS")[k]
    cols<-list(sex_col,eth_col,iris_col)[[k]]
    level_order<-list(c("Male","Female"),c("Unknown","Hispanics","Caucasian","Black","Asian"),c("Unknown","IS","IR"))[[k]]
    tab<-table(factor(d[[variable]],levels=level_order)); y0<-.03
    for(j in seq_along(tab)){
      y1<-y0+tab[j]/nrow(d)*.66
      rect(x-.05,y0,x+.05,y1,col=cols[names(tab)[j]],border="grey35",lwd=.4)
      if(tab[j]/nrow(d)>.12)text(x,(y0+y1)/2,tab[j],cex=5/7,col="grey20")
      y0<-y1
    }
    text(x,.725,variable,cex=5/7,srt=90,adj=0)
  }
},file.path(fig,"01_cohort_multitrack"),180,180)

## 2. 可独立复用的年龄与BMI摘要，各自一个文件 ====
for(v in c("age","BMI")) {
  lim<-if(v=="age") c(20,80) else c(16,44)
  col<-if(v=="age") age_col else bmi_col
  save_base({
  par(ps=11)
    par(mar=c(1,3,1,1),mgp=c(1.6,.35,0),tcl=-.15)
    plot(1,d[[v]][1],type="n",xlim=c(.5,1.5),ylim=lim,xlab="",ylab=if(v=="age")"Age (years)" else "BMI (kg/m2)",xaxt="n",yaxt="n",bty="l",cex.axis=5/7,cex.lab=6/7)
    axis(2,at=if(v=="age")seq(20,80,10)else seq(20,40,5),cex.axis=5/7,lwd=.5,lwd.ticks=.5)
    b<-round(d[[v]]/(diff(lim)/55)); off<-ave(seq_along(b),b,FUN=function(z)(seq_along(z)-(length(z)+1)/2)*.028)
    raster_points(1+off,d[[v]],col=col,pch=16,cex=.4)
    boxplot(d[[v]],at=1,add=TRUE,axes=FALSE,outline=FALSE,boxwex=.4,col=NA,lwd=.5)
  },file.path(fig,paste0("01_",v,"_summary")),45,85)
}
write.csv(data.frame(variable=c("age","BMI"),n=nrow(d),mean=c(mean(d$age),mean(d$BMI)),sd=c(sd(d$age),sd(d$BMI))),file.path(result,"01_cohort_summaries.csv"),row.names=FALSE)

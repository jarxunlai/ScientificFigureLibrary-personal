# 表型九类指标渐变柱与分组折线 / Phenotype profiles with gradient bars and grouped lines
# 科学问题：不同表型在九类标准化特征上的偏离模式是否一致？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-phenotype-bar-line")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：放大正负值相对零参考环的起伏，恢复渐变柱、误差和四表型折线的独立版本；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 九组特征：渐变条形与四表型折线。每个版本独立文件。
## ci_half是直接生成的示例95%CI半宽；星号来自输入，不据此声称统计发现。
#############################################################################


source(file.path(module_root,"code/export.R"))
bars<-read.csv(file.path(module_root,"data/phenotype_bars.csv"));lines_data<-read.csv(file.path(module_root,"data/phenotype_lines.csv"))
fig<-file.path(entry,"figure");groups<-unique(bars$group)
counts<-as.numeric(table(factor(bars$group,levels=groups)))
step<-360/(nrow(bars)+2*length(groups));bars$angle<-NA_real_
start<-90
for(k in seq_along(groups)){ix<-which(bars$group==groups[k]);bars$angle[ix]<-start-(seq_along(ix)-.5)*step;start<-start-(length(ix)+2)*step}
lines_data$angle<-bars$angle[match(lines_data$feature,bars$feature)]
gradient<-circlize::colorRamp2(c(-1.5,0,1.5),c("#4575B4","#FFFFBF","#F34324"))
line_colors<-c(PT1="#D45E00",PT2="#E59F00",PT3="#CA78A2",PT4="#0070B3")
rzero<-.68;scale<-.17
write.csv(bars,file.path(entry,"result/08_phenotype_bars.csv"),row.names=FALSE)
write.csv(lines_data,file.path(entry,"result/08_phenotype_lines.csv"),row.names=FALSE)

## 1. 共用角度映射，图层代码内联；参考零值是半径0.70的圆而非圆心 ====
for(kind in c("gradient_bars","grouped_lines")){
  save_base({
  par(ps=11)
    canvas(c(-1.28,1.38),c(-1.28,1.28))
    points_xy<-matrix(numeric(0),ncol=2);points_col<-character()
    sector(0,360,0,.425,if(kind=="gradient_bars")"#F6DECC" else "#F9F9F9",NA)
    for(k in seq_along(groups)){
      ix<-which(bars$group==groups[k]);a0<-max(bars$angle[ix])+step/2;a1<-min(bars$angle[ix])-step/2
      arc(a1,a0,.43,"black",.6);radial_segment(mean(c(a0,a1)),.415,.43,lwd=.55)
      # 分组短线与内圈名称；九个组名字体和表型标题分开控制。
      pos<-xy(mean(c(a0,a1)),.325)
      text(pos[1,1],pos[1,2],gsub(" ","\n",groups[k]),cex=5/7)
      for(v in seq(-1.5,1.5,.5)){
        next_a<-a1-2*step
        arc(next_a,a1,rzero+v*scale,"grey60",.45)
      }
      arc(a1,a0,rzero,"grey30",.7)
      if(kind=="grouped_lines"){
        for(pt in names(line_colors)){
          z<-lines_data[lines_data$group==groups[k]&lines_data$phenotype==pt,]
          z<-z[match(bars$feature[ix],z$feature),]
          p<-xy(z$angle,rzero+z$value*scale);lines(p,col=line_colors[pt],lwd=.7)
          points_xy<-rbind(points_xy,p);points_col<-c(points_col,rep(line_colors[pt],nrow(p)))
        }
      }
      for(i in ix){
        if(kind=="gradient_bars"){
          rr<-rzero+bars$mean[i]*scale;e<-bars$ci_half[i]*scale
          sector(bars$angle[i]-.33*step,bars$angle[i]+.33*step,rzero,rr,gradient(bars$mean[i]),NA)
          radial_segment(bars$angle[i],rr-e,rr+e,"grey50",.5)
          arc(bars$angle[i]-.7,bars$angle[i]+.7,rr+e,"grey50",.5)
          labr<-max(rzero,rr+e)+.025
        } else {
          labr<-max(rzero+lines_data$value[lines_data$feature==bars$feature[i]]*scale)+.025
        }
        radial_text(bars$angle[i],labr,paste0(bars$feature[i],bars$significance[i]),"grey40",5/7)
      }
    }
    if(kind=="grouped_lines")raster_points(points_xy,pch=16,col=points_col,cex=.29)
    for(v in seq(-1.5,1.5,.5)){p<-xy(90,rzero+v*scale);text(p[1,1],p[1,2],if(v>0)paste0("+",v)else as.character(v),cex=5/7)}
    if(kind=="gradient_bars"){
      text(0,0,"PT 2\npsychosomatic\nexample profile",cex=6/7,font=2)
      color_key(gradient,c(-1.5,1.5),1.15,.76,height=.35,ticks=c(-1.5,0,1.5),title="z score")
    } else {
      legend(.72,1.20,legend=c("PT 1: avoidant","PT 2: psychosomatic","PT 3: somatic","PT 4: distress"),col=line_colors,lty=1,pch=16,cex=5/7,bty="n",y.intersp=.9)
    }
  },file.path(fig,paste0("08_tinnitus_",kind)),180,180)
}

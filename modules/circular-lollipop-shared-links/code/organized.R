# 疾病蛋白选择得分与共享连线 / Disease protein selection scores and shared links
# 科学问题：哪些蛋白同时出现在多个疾病的预测特征集合中？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-lollipop-shared-links")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复每个疾病内得分递减的棒棒糖曲线及同名蛋白共享连线；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 跨疾病蛋白选择得分：12疾病等宽、组内蛋白均匀分布。
## 内圈线只连输入表的同名蛋白；score在0~1，不作因果解释。
#############################################################################


source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/protein_scores.csv"));fig<-file.path(entry,"figure")
groups<-unique(d$disease);d<-d[order(match(d$disease,groups),-d$score),];span<-340/length(groups);gap<-2
d$angle<-NA_real_
for(k in seq_along(groups)){
  ix<-which(d$disease==groups[k]);d$angle[ix]<-90-(k-1)*span-gap/2-(seq_along(ix)-.5)*(span-gap)/length(ix)
}
d$color<-ifelse(d$shared,"#3690BC","#F5A991")
links<-do.call(rbind,lapply(unique(d$protein[d$shared]),function(g){ix<-which(d$protein==g);pairs<-combn(ix,2);data.frame(protein=g,from=pairs[1,],to=pairs[2,])}))
write.csv(d,file.path(entry,"result/03_protein_scores.csv"),row.names=FALSE)
write.csv(links,file.path(entry,"result/03_shared_protein_links.csv"),row.names=FALSE)

## 1. 标准化棒棒糖、疾病环和共享蛋白内连线 ====
save_base({
  par(ps=11)
  canvas(c(-1.18,1.18),c(-1.18,1.18))
  for(i in seq_len(nrow(links)))curve_link(d$angle[links$from[i]],.43,d$angle[links$to[i]],.43,.03,adjustcolor("#4E6A74",alpha.f=.45),.4)
  for(k in seq_along(groups)){
    a0<-90-(k-1)*span-gap/2;a1<-90-k*span+gap/2
    sector(a1,a0,.84,.91,"#EEEEEE","grey40",.45)
    arc_text((a0+a1)/2,.873,groups[k],cex=5/7)
    sector(a1,a0,.43,.825,"white","grey65",.4)
    for(t in c(.2,.4,.6,.8))arc(a1,a0,.43+t*.395,"grey85",.35)
    ix<-which(d$disease==groups[k])
    for(i in ix){rr<-.43+d$score[i]*.395;radial_segment(d$angle[i],.43,rr,d$color[i],.65);p<-xy(d$angle[i],rr);points(p,pch=16,col=d$color[i],cex=.42)
      radial_segment(d$angle[i],.915,.944,"grey60",.35);radial_text(d$angle[i],.95,d$protein[i],d$color[i],5/7)}
  }
  a0<-90-gap/2
  for(t in seq(0,1,.2)){p<-xy(a0,.43+t*.395);text(p[1,1]-.025,p[1,2],format(t,nsmall=1),cex=5/7,adj=1)}
  text(-.09,.62,"Selection score",srt=90,cex=6/7)
  legend("bottom",inset=.015,legend=c("Shared protein","Disease-specific protein"),col=c("#3690BC","#F5A991"),pch=16,bty="n",horiz=TRUE,cex=5/7)
},file.path(fig,"03_protein_disease_lollipop"),180,180)

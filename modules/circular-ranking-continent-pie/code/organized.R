# 国家排名外圈与洲别贡献饼图 / Country ranking bars with a continent contribution pie
# 科学问题：国家贡献排名与洲别总体占比如何保持数量一致？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-ranking-continent-pie")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- .48 # 保留已审阅模板的数值径向基线；不是文件路径别名。
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复国家排名标签、径向刻度与中心洲别饼图的视觉比例，四主题仍独立；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 国家排名与洲别占比：饼图严格聚合同一次绘制的国家数值。
## 分别导出教程紫色版、来源橙蓝绿三主题；不拼成多panel文件。
#############################################################################


source(file.path(module_root,"code/export.R"))
raw<-read.csv(file.path(module_root,"data/country_contributions.csv"),na.strings=character());fig<-file.path(entry,"figure")
regions<-c("AF","AS","EU","NA","OA","SA")
themes<-c("tutorial_purple","source_orange","source_blue","source_green")
metrics<-c("Value","Value","Urban","Rural")
palettes<-list(c("#E0D3EC","#BEA4DB","#9D7CC5","#8160AD","#68429A","#4D267F"),
 c("#F5C9BF","#F0AD9B","#EB8963","#DF6C2E","#C65F28","#AB5121"),
 c("#D3E8F2","#A6D1E4","#79B8D1","#4C9EBA","#297EA0","#16647F"),
 c("#E4F1D6","#CAE6B1","#ABD386","#83BE5C","#5AA443","#2F8734"))
bar_cols<-c("#7241A2","#EF7331","#268FAE","#58AD43")

## 1. 每个国家一个角度单元，不因中心饼图而改变国家排名 ====
for(k in seq_along(themes)){
  d<-raw;d$PlotValue<-d[[metrics[k]]];d<-d[order(d$PlotValue,decreasing=TRUE),]
  pie<-aggregate(PlotValue~Region,d,sum);pie<-pie[match(regions,pie$Region),]
  pie$percentage<-pie$PlotValue/sum(pie$PlotValue)*100
  write.csv(d,file.path(entry,"result",paste0("10_",themes[k],"_countries.csv")),row.names=FALSE)
  write.csv(pie,file.path(entry,"result",paste0("10_",themes[k],"_continents.csv")),row.names=FALSE)
  save_base({
  par(ps=11)
    canvas(c(-1.12,1.16),c(-1.10,1.23))
    step<-350/nrow(d);max_value<-ceiling(max(d$PlotValue))
    for(i in seq_len(nrow(d))){
      a<-90-(i-.5)*step;rr<-base+d$PlotValue[i]/max_value*.56
      sector(a-.44*step,a+.44*step,base,rr,bar_cols[k],"white",.25)
      radial_text(a,rr+.015,d$Country[i],cex=5/7)
    }
    radial_segment(90,base,base+.56,lwd=.6)
    for(v in 0:max_value){y<-base+v/max_value*.56;segments(-.012,y,.012,y,lwd=.5);text(-.023,y,v,cex=5/7,adj=1)}
    text(.035,1.08,paste0(metrics[k]," contribution"),srt=90,adj=0,cex=5/7)
    # 中心扇形与外圈共用d的聚合，不用独立编造百分比。
    st<-90;palette<-setNames(palettes[[k]],regions)
    for(j in seq_len(nrow(pie))){
      en<-st-pie$percentage[j]/100*360
      sector(en,st,0,.405,palette[pie$Region[j]],"white",.45)
      mid<-(st+en)/2
      if(pie$percentage[j]>=5){
        p<-xy(mid,.245);text(p[1,1],p[1,2],sprintf("%s\n%.1f%%",pie$Region[j],pie$percentage[j]),cex=5/7)
      } else {
        p<-xy(mid,.42);p2<-xy(mid,.455);segments(p[1,1],p[1,2],p2[1,1],p2[1,2],lwd=.4)
        radial_text(mid,.46,sprintf("%s %.1f%%",pie$Region[j],pie$percentage[j]),cex=5/7)
      }
      st<-en
    }
  },file.path(fig,paste0("10_country_continent_",themes[k])),180,180)
}

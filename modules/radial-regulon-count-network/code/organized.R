# 调控基因集放射网络与外圈数量 / Radial regulon network with outer counts
# 科学问题：不同功能类别下的基因集规模及层级归属如何同时展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "radial-regulon-count-network")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：缩短过度突出的外圈条形、重排类别多行名称，恢复中心分支与外圈数量的层级；输入/派生值均为合成示例，不冒充来源结果。
#############################################################################
## 层级网络与外圈数量：相同坐标绑定边、叶节点、标签和条形。
## 输入是synthetic分组表；边是明确层级关系，不是重建的调控作用。
#############################################################################


source(file.path(module_root,"code/export.R"))
d<-read.csv(file.path(module_root,"data/gene_sets.csv"));fig<-file.path(entry,"figure")
categories<-unique(d$category)
cols<-setNames(c("#F76502","#F41237","#026026","#2F0D06","#0372E9","#F00A0B","#2E3426","#DEA34A","#AA1F19","#5728A0"),categories)
d<-d[order(match(d$category,categories),-d$number),]
step<-360/nrow(d);d$angle<-90-(seq_len(nrow(d))-.5)*step
centers<-aggregate(angle~category,d,mean);centers<-centers[match(categories,centers$category),]
write.csv(d,file.path(entry,"result/05_gene_sets_layout.csv"),row.names=FALSE)

## 1. 完整层级绘图可见；条形版和气泡版各自独立导出 ====
for(style in c("bars","bubbles")){
  save_base({
  par(ps=11)
    canvas(c(-1.25,1.25),c(-1.25,1.25))
    for(k in seq_along(categories)){
      g<-categories[k];a<-centers$angle[k];node<-xy(a,.37)
      radial_segment(a,0,.37,adjustcolor(cols[g],.65),.6)
      points(node,pch=16,cex=.75,col=cols[g])
      ix<-which(d$category==g)
      for(i in ix) {
        curve_link(a,.37,d$angle[i],.83,.48,adjustcolor(cols[g],.34),.45)
        if(style=="bars")sector(d$angle[i]-.37*step,d$angle[i]+.37*step,.83,.83+d$number[i]/max(d$number)*.13,adjustcolor(cols[g],.75),NA)
        if(style=="bubbles"){p<-xy(d$angle[i],.87);points(p,pch=21,bg=adjustcolor(cols[g],.65),col=cols[g],cex=sqrt(d$number[i]/max(d$number))*1.6)}
        rr<-if(style=="bars") .845+d$number[i]/max(d$number)*.13 else .92
        radial_text(d$angle[i],rr,d$set[i],cols[g],5/7)
      }
      # 短行类别名放在类别节点内侧，避免把生物学参数藏进helper。
      label<-paste(strwrap(g,width=16),collapse="\n")
      pos<-xy(a,.37);pos2<-xy(a,.47)
      text(pos2[1,1],pos2[1,2],label,cex=.75,col=cols[g],font=2)
    }
    points(0,0,pch=16,cex=.5,col="grey40")
    if(style=="bars"){
      text(0,-1.20,"Outer bars: number of genes",cex=5/7)
    } else {
      legend("bottom",inset=.01,legend=c("10 genes","40 genes","70 genes"),pt.cex=sqrt(c(10,40,70)/max(d$number))*1.6,pch=21,pt.bg="grey70",bty="n",horiz=TRUE,cex=5/7)
    }
  },file.path(fig,paste0("05_regulon_network_",style)),180,180)
}

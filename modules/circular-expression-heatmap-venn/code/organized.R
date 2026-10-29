# 环状分组表达热图与集合交集 / Circular grouped expression heatmap with Venn inset
# 科学问题：不同形态发生相关基因集合在组织区域中的表达及集合重叠如何？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "circular-expression-heatmap-venn")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复基因名、表达和q轨道的可读比例及中心交集标记，保留真实示例成员计数；输入/派生值均为合成示例，不冒充来源结果。
# 环状分组表达热图 + qval 菱形 + 中心集合交集。
# 输入：独立模拟表达与集合标签；组内按 qval 排序，表达按组/区域缩放到 0~0.5。
# 中心为固定圆形集合示意，面积不编码数量；全部 7 个区域计数来自同一输入成员集合。
library(circlize)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)


source(file.path(module_root,"code/export.R"))

## 1. 数据、排序、表达缩放和可追溯集合计数 ====
d<-read.csv(file.path(module_root,"data/circular_expression.csv"),check.names=FALSE)
group_order<-c("Divergence","Curl","Acceleration")
d$Group<-factor(d$Group,levels=group_order)
d<-d[order(d$Group,d$qval),]
regions<-c("Left atrium","Right atrium","Outflow tract","Right ventricle","Left ventricle")
for(g in group_order) for(v in regions) {
  idx<-which(d$Group==g)
  d[idx,v]<-scales::rescale(d[idx,v],to=c(0,.5))
}
m<-as.matrix(d[,regions]); rownames(m)<-d$Gene
sets<-split(d$Gene,d$Group)
all_genes<-sort(unique(d$Gene))
membership<-sapply(sets,function(s) all_genes %in% s)
codes<-apply(membership,1,paste0,collapse="")
counts<-table(factor(codes,levels=c("TRUEFALSEFALSE","FALSETRUEFALSE","FALSEFALSETRUE","TRUETRUEFALSE","TRUEFALSETRUE","FALSETRUETRUE","TRUETRUETRUE")))
write.csv(data.frame(Gene=all_genes,membership),file.path(draft,"result/06_set_membership.csv"),row.names=FALSE)
write.csv(data.frame(Region=names(counts),Count=as.integer(counts)),file.path(draft,"result/06_venn_counts.csv"),row.names=FALSE)
write.csv(d,file.path(draft,"result/06_normalized_expression.csv"),row.names=FALSE)
col_fun<-colorRamp2(seq(0,.5,length.out=5),rev(brewer.pal(5,"PiYG")))
q_col<-colorRamp2(c(0,.0025,.005,.0075,.01),c("#fff1e7","#f99c70","#dd4139","#722043","black"))
g_col<-c(Divergence="#4DAC26",Curl="#F7F7F7",Acceleration="#D01C8B")

## 2. 环形轨道：外侧基因名、表达、qval、内侧组名 ====
save_draw({
  par(mar=c(0,0,0,0),family="Arial")
  circos.clear()
  circos.par(start.degree=60,gap.after=c(5,5,30),track.margin=c(0,.008),cell.padding=c(0,0,0,0),
             canvas.xlim=c(-1.17,1.33),canvas.ylim=c(-1.15,1.15),points.overflow.warning=FALSE)
  circos.heatmap(m,split=d$Group,cluster=FALSE,col=col_fun,track.height=.24,
                 bg.border="black",bg.lwd=.6,cell.border="white",cell.lwd=.35,
                 rownames.side="outside",rownames.cex=.66,rownames.font=1)
  # 列名位于开口处；从内外顺序反向匹配环形矩阵。
  circos.track(track.index=get.current.track.index(),bg.border=NA,panel.fun=function(x,y){
    if(CELL_META$sector.numeric.index==3) {
      ypos<-seq(CELL_META$cell.ylim[1],CELL_META$cell.ylim[2],length.out=6)
      ypos<-(head(ypos,-1)+tail(ypos,-1))/2
      circos.text(rep(CELL_META$xlim[2],5)+convert_x(1.5,"mm"),ypos,rev(regions),
                  facing="inside",adj=c(0,.5),cex=.64)
    }
  })
  circos.track(ylim=c(0,1),track.height=.045,bg.border=NA,panel.fun=function(x,y){
    block<-d[d$Group==CELL_META$sector.index,]
    circos.points(seq_len(nrow(block))-.5,rep(.5,nrow(block)),pch=18,cex=.55,col=q_col(block$qval))
  })
  circos.track(ylim=c(0,1),track.height=.075,bg.col=adjustcolor(g_col,alpha.f=.3),bg.border="#666666",bg.lwd=.5,
    panel.fun=function(x,y){
      circos.text(CELL_META$xcenter,.5,CELL_META$sector.index,facing="bending.inside",niceFacing=TRUE,cex=.82)
    })
  # 中心固定圆形示意，精确计数而非假设圆面积按比例。
  pushViewport(viewport(x=.477,y=.53,width=.29,height=.27))
  grid.circle(.32,.58,r=.23,gp=gpar(fill=adjustcolor(g_col[1],alpha.f=.28),col="#555555",lty=2,lwd=.5))
  grid.circle(.55,.58,r=.23,gp=gpar(fill=adjustcolor("#cccccc",alpha.f=.22),col="#555555",lty=2,lwd=.5))
  grid.circle(.435,.37,r=.23,gp=gpar(fill=adjustcolor(g_col[3],alpha.f=.28),col="#555555",lty=2,lwd=.5))
  pos<-rbind(c(.19,.63),c(.68,.63),c(.435,.23),c(.435,.67),c(.30,.43),c(.57,.43),c(.435,.50))
  grid.text(as.integer(counts),x=pos[,1],y=pos[,2],gp=gpar(fontsize=8,fontfamily="Arial"))
  grid.text("Morphogenic\ngenes",x=.435,y=-.06,gp=gpar(fontsize=10,fontface="bold",fontfamily="Arial"))
  popViewport()
  lg<-Legend(title="Expression\nfraction",col_fun=col_fun,at=c(0,.25,.5),grid_height=unit(2.4,"mm"),
             legend_height=unit(23,"mm"),labels_gp=gpar(fontsize=8),title_gp=gpar(fontsize=8))
  draw(lg,x=unit(.965,"npc"),y=unit(.92,"npc"),just=c("right","top"))
  lg<-Legend(title="q value",col_fun=q_col,at=c(0,.005,.01),legend_height=unit(23,"mm"),
             labels_gp=gpar(fontsize=8),title_gp=gpar(fontsize=8))
  draw(lg,x=unit(.965,"npc"),y=unit(.14,"npc"),just=c("right","bottom"))
  circos.clear()
},file.path(draft,"figure/06_circular_venn"),180,176)

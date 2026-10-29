# 环形多维注释UMAP / UMAP with circular multi-dimensional annotations
# 科学问题：细胞嵌入、类型组成、时间点和条件信息如何在同一图中对应？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "umap-circular-multitrack")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 环形多轨 UMAP：来源教程中的中心非高斯形态、外圈 log10 刻度、三轨与图例。
# 当前 data/cell_embedding.csv 为合成示例；本轮形态/数量/元数据同样为合成。
# 不执行 Seurat / UMAP / Harmony；不读取来源 PNG 作为点层，不复制论文坐标。
# 所有真正绘图输入、环比例、数量都输出 result/ 可审查。
Sys.setlocale("LC_CTYPE","English_United States.utf8")



source(file.path(module_root,"code/export.R"))
library(circlize)
fig <- file.path(base,"figure"); res <- file.path(base,"result")
raw <- read.csv(file.path(module_root,"data/cell_embedding.csv"))
set.seed(260119)
types <- c("Amacrine_cell","Astrocyte","B_cell","Bipolar_cell","Endothelial","MG",
           "Microglia","Monocyte_macrophage","Neutrophil","RGC","RPE","T_cell")
cols <- setNames(c("#087ca5","#4e4d9a","#00b0cf","#00a3a3","#d70135","#ef6660",
                  "#009e58","#c6df71","#8b1697","#ffca00","#dfa924","#f65ea5"),types)
time_cols <- c(Control="#f15c4d","12h"="#ffd966","1d"="#00a88b","2d"="#17afc9","4d"="#494492","1w"="#80bb48","2w"="#ff840b")
condition_cols <- c(Control="#2c7fb8",Injury="#f7950b")
# 手工设计的示例轮廓：非来源图像像素/真实数据的数值提取。
shapes <- list(
  Amacrine_cell=rbind(c(-.68,-.30),c(-.66,-.18),c(-.60,-.14),c(-.52,-.23),c(-.48,-.38),c(-.55,-.43),c(-.65,-.39)),
  Astrocyte=rbind(c(-.41,-.17),c(-.39,-.03),c(-.43,.08),c(-.32,.13),c(-.21,.10),c(-.11,.15),c(.01,.06),c(.10,.06),c(.13,-.02),c(.06,-.05),c(.15,-.15),c(.07,-.21),c(-.08,-.19),c(-.18,-.25),c(-.30,-.23)),
  B_cell=rbind(c(.30,-.58),c(.40,-.50),c(.44,-.49),c(.48,-.53),c(.37,-.55),c(.31,-.60)),
  Bipolar_cell=rbind(c(-.71,.08),c(-.69,.13),c(-.67,.10),c(-.68,.04)),
  Endothelial=rbind(c(.39,-.16),c(.36,-.10),c(.36,-.01),c(.39,.02),c(.40,-.06),c(.43,-.13)),
  MG=rbind(c(-.28,.25),c(-.33,.38),c(-.32,.53),c(-.22,.68),c(-.16,.79),c(.02,.85),c(.14,.82),c(.22,.73),c(.25,.66),c(.36,.62),c(.38,.53),c(.30,.46),c(.26,.32),c(.12,.30),c(.03,.35),c(-.11,.28)),
  Microglia=rbind(c(-.36,-.65),c(-.24,-.57),c(-.11,-.62),c(.01,-.55),c(.05,-.42),c(.12,-.39),c(.16,-.53),c(.22,-.62),c(.20,-.68),c(.28,-.77),c(.21,-.84),c(.06,-.82),c(-.05,-.86),c(-.18,-.80),c(-.28,-.78)),
  Monocyte_macrophage=rbind(c(-.52,-.68),c(-.51,-.59),c(-.40,-.55),c(-.37,-.47),c(-.30,-.50),c(-.24,-.45),c(-.21,-.57),c(-.26,-.64),c(-.22,-.73),c(-.33,-.74),c(-.38,-.67)),
  Neutrophil=rbind(c(-.53,-.78),c(-.50,-.73),c(-.47,-.77),c(-.46,-.85),c(-.51,-.85)),
  RGC=rbind(c(-.76,-.34),c(-.72,-.28),c(-.69,-.27),c(-.66,-.36),c(-.61,-.39),c(-.60,-.44),c(-.68,-.43),c(-.73,-.40)),
  RPE=rbind(c(.57,.26),c(.66,.28),c(.73,.22),c(.77,.12),c(.78,.01),c(.82,-.07),c(.83,-.20),c(.80,-.32),c(.75,-.29),c(.72,-.20),c(.69,-.13),c(.68,-.03),c(.64,.08),c(.57,.12),c(.54,.20)),
  T_cell=rbind(c(-.67,-.64),c(-.63,-.65),c(-.57,-.70),c(-.54,-.73),c(-.59,-.74),c(-.65,-.70)))
counts <- setNames(c(1800,10000,500,350,480,22500,12000,3400,600,900,5700,650),types)
clouds <- list()

## 1. 使用现有合成 CSV 提供秩扰动，恢复非高斯轮廓、空洞与细胞数量差异 ====
for(k in seq_along(types)) {
  ct <- types[k]; poly <- shapes[[ct]]; n <- counts[ct]
  original <- raw[raw$Celltype==ct,]
  # 均匀候选经原 CSV 秩扰动；不是原始单细胞嵌入重计算。
  u0 <- (rank(original$UMAP1)-.5)/nrow(original)
  v0 <- (rank(original$UMAP2)-.5)/nrow(original)
  nn <- n*8
  u <- (runif(nn)+rep(u0,length.out=nn))%%1
  v <- (runif(nn)+rep(v0,length.out=nn))%%1
  tt <- seq(0,1,length.out=nrow(poly)+1)
  px <- spline(tt,c(poly[,1],poly[1,1]),method="periodic",n=250)$y
  py <- spline(tt,c(poly[,2],poly[1,2]),method="periodic",n=250)$y
  x <- min(px)+u*diff(range(px));y <- min(py)+v*diff(range(py))
  keep <- sp::point.in.polygon(x,y,px,py)>0
  if(ct=="MG") {
    # 两处长条空洞形成教程中不规则的上方大片细胞云。
    keep <- keep & !(((x+.02)/.19)^2+((y-.44)/.033)^2<1 & runif(nn)>.11) &
                   !(((x-.18)/.11)^2+((y-.57)/.018)^2<1 & runif(nn)>.18)
  }
  if(ct=="Astrocyte") keep <- keep & !(((x+.18)/.11)^2+((y+.045)/.018)^2<1 & runif(nn)>.10)
  take <- which(keep)[seq_len(n)]
  d <- data.frame(Celltype=ct,x=x[take]+rnorm(n,sd=.005),y=y[take]+rnorm(n,sd=.005))
  # 混合由同一份合成元数据生成，并输出；Control在MG示例富集，不声称真实比例。
  time_prob <- if(ct=="MG") c(.78,.06,.03,.03,.03,.03,.04) else c(.08,.07,.15,.18,.18,.17,.17)
  d$Timepoint <- sample(names(time_cols),n,replace=TRUE,prob=time_prob)
  injury_p <- if(ct %in% c("MG","RPE","Microglia")) .83 else .62
  d$Condition <- ifelse(runif(n)<injury_p,"Injury","Control")
  d$Cluster <- sprintf("Cluster%02d",k)
  clouds[[ct]] <- d
}
df <- bind_rows(clouds)
write.csv(df,file.path(res,"19_visual_example_embedding.csv"),row.names=FALSE)
sector_counts <- data.frame(Celltype=types,n=as.integer(counts),log10_n=log10(counts))
write.csv(sector_counts,file.path(res,"19_visual_example_sector_counts.csv"),row.names=FALSE)
ring <- list()
for(ct in types) for(track in c("Timepoint","Condition")) {
  tab <- prop.table(table(df[df$Celltype==ct,track]))
  ring[[length(ring)+1]] <- data.frame(Celltype=ct,track=track,group=names(tab),proportion=as.numeric(tab))
}
write.csv(bind_rows(ring),file.path(res,"19_visual_example_ring_proportions.csv"),row.names=FALSE)

## 2. 用 circlize 恢复弯曲标签/数量刻度/开口轨道名称；四种中心着色保持同一几何 ====
for(view in c("Celltype","Timepoint","Condition","Cluster")) {
  point_pal <- switch(view,Celltype=cols,Timepoint=time_cols,Condition=condition_cols,
                      Cluster=setNames(unname(cols),sprintf("Cluster%02d",1:12)))
  for(device in c("png","pdf")) {
    stem <- file.path(fig,paste0("19_circular_umap_",view))
    if(device=="png") ragg::agg_png(paste0(stem,".png"),width=240,height=185,units="mm",res=300,background="white")
    if(device=="pdf") cairo_pdf(paste0(stem,".pdf"),width=240/25.4,height=185/25.4,family="Arial",bg="white")
    par(mar=c(0,0,0,0),family="Arial",ps=10,xpd=NA)
    circos.clear()
    circos.par(cell.padding=c(0,0,0,0),track.margin=c(.006,0),
      gap.degree=c(rep(2,11),12),start.degree=0,
      canvas.xlim=c(-1.20,1.91),canvas.ylim=c(-1.20,1.20),points.overflow.warning=FALSE)
    circos.initialize(factor(types,levels=types),xlim=cbind(rep(0,12),log10(counts)))
    circos.trackPlotRegion(ylim=c(0,1),track.height=.012,bg.border=NA,panel.fun=function(x,y) {
      ct <- CELL_META$sector.index
      circos.text(CELL_META$xcenter,1+mm_y(4.1),ct,cex=1.02,facing="bending.inside",niceFacing=TRUE)
      ticks <- 0:floor(log10(counts[ct]))
      circos.axis(h="top",major.at=ticks,labels=ticks,labels.cex=.46,minor.ticks=1,
        major.tick.length=mm_y(.7),labels.niceFacing=TRUE,lwd=.4)
      circos.rect(0,0,log10(counts[ct]),1,col=cols[ct],border=NA)
    })
    for(track in c("Timepoint","Condition")) {
      pal <- if(track=="Timepoint") time_cols else condition_cols
      circos.trackPlotRegion(ylim=c(0,1),track.height=.010,bg.border=NA,panel.fun=function(x,y) {
        ct <- CELL_META$sector.index
        tab <- table(factor(df[df$Celltype==ct,track],levels=names(pal)))
        ends <- cumsum(tab)/sum(tab)*log10(counts[ct]);starts <- c(0,head(ends,-1))
        circos.rect(starts,0,ends,1,col=pal,border=NA)
      })
    }
    # 中心数据坐标本身是此脚本生成的规范化示例，不再被包自动放大为等间隔圆点。
    points(df$x,df$y,pch=16,cex=.17,col=scales::alpha(point_pal[df[[view]]],.63))
    text(c(1.032,.997,.964),c(.075,.075,.075),c("Celltype","Timepoint","Condition"),srt=-90,cex=.65)
    legend_sets <- list(Celltype=cols,Timepoint=time_cols,Condition=condition_cols)
    if(view=="Cluster") names(legend_sets$Celltype) <- paste(sprintf("C%02d",1:12),types)
    y <- .75
    for(lg in names(legend_sets)) {
      pp <- legend_sets[[lg]]
      text(1.20,y,lg,adj=c(0,.5),cex=.98)
      for(i in seq_along(pp)) {
        yy <- y-i*.049
        rect(1.20,yy-.018,1.247,yy+.018,col=pp[i],border=NA)
        text(1.26,yy,names(pp)[i],adj=c(0,.5),cex=.82)
      }
      y <- y-(length(pp)+1.7)*.049
    }
    circos.clear(); dev.off()
  }
}
cat("VISUAL_REPAIR:",nrow(df),"synthetic points; four views; 3 rings with count axes.\n")

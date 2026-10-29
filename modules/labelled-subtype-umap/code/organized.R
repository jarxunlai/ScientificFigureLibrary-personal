# T细胞亚群标注UMAP / Labelled T-cell-subtype UMAP
# 科学问题：T细胞亚群在已有嵌入坐标上如何分布并保持标签可读？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "labelled-subtype-umap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 带标签 UMAP：恢复来源两幅图（全细胞群 + T细胞亚群），不缺省第一幅。
# data/tcell_embedding.csv 是合成示例，本轮仍只生成二维拟态示例、标签及配色。
# 无Seurat对象、无真实重聚类/Harmony/UMAP；不把图形示例当作Cell论文结果。
# 轮廓是人工设计，不读取来源图片作为散点层；绘图输入完整保存在 result/。
Sys.setlocale("LC_CTYPE","English_United States.utf8")



source(file.path(module_root,"code/export.R"))
fig <- file.path(base,"figure");res <- file.path(base,"result")
raw <- read.csv(file.path(module_root,"data/tcell_embedding.csv"))
set.seed(260123)

## 1. T细胞：单个连续轮廓与软边界分区，而不是十个分离的椭圆团 ====
cols <- c("CD8 Chemokine producing"="#7dff56","CD8 IFN-responsive"="#c9fc36",
  "CD4 Tissue resident"="#4453c4","CD8 Undefined"="#cccccc","CD8 Proliferative"="#f1ca3a",
  "CD8 Tissue resident"="#bf2101","CD8 Early activated"="#ea4f0d","CD8 Stem-like"="#2bf6a7",
  "CD8 NK-like"="#fe922a","CD4 T regulatory"="#4490fe")
centers <- data.frame(Celltype=names(cols),
  x=c(1.0,.6,3.8,.9,-2.6,-4.6,-.5,-.3,-6.2,4.9),
  y=c(-1.1,1.3,2.25,-2.8,-3.7,-.3,-.25,3.05,1.4,-.55),
  sx=c(1.5,1.7,1.4,1.5,1.1,1.7,1.3,1.8,.9,1.3),
  sy=c(.7,.8,1.2,1.4,.9,1.6,.8,1,.8,1.0))
poly <- rbind(c(-8.1,1.65),c(-7.4,2.25),c(-6.3,2.42),c(-5.3,2.75),c(-3.8,3.50),
  c(-1.7,4.10),c(.3,4.25),c(2.5,3.95),c(4.5,3.4),c(5.55,2.5),c(5.7,1.5),
  c(5.0,.55),c(6.2,.1),c(6.65,-.8),c(6.1,-1.75),c(4.7,-2.4),c(3.5,-2.55),
  c(2.8,-3.65),c(1.3,-4.35),c(-.5,-4.2),c(-2.6,-4.8),c(-3.5,-4.4),
  c(-4.2,-3.25),c(-5.8,-2.05),c(-6.6,-.95),c(-7.8,.05))
# splined outline retains a connected cell cloud; not 10 Gaussian islands.
tt <- seq(0,1,length.out=nrow(poly)+1)
outline <- list(x=spline(tt,c(poly[,1],poly[1,1]),method="periodic",n=350)$y,
                y=spline(tt,c(poly[,2],poly[1,2]),method="periodic",n=350)$y)
u <- (rank(raw$UMAP1)-.5)/nrow(raw);v <- (rank(raw$UMAP2)-.5)/nrow(raw)
nn <- 500000
x <- -8.5+15.5*((runif(nn)+rep(u,length.out=nn))%%1)
y <- -5+9.7*((runif(nn)+rep(v,length.out=nn))%%1)
inside <- sp::point.in.polygon(x,y,outline$x,outline$y)>0
x <- x[inside][1:170000];y <- y[inside][1:170000]
# 空间相邻标签混合：来源的连续状态区，而不是一套新生物学分类。
score <- sapply(1:nrow(centers),function(k)
  -((x-centers$x[k])/centers$sx[k])^2-((y-centers$y[k])/centers$sy[k])^2+rnorm(length(x),sd=.55))
ct <- centers$Celltype[max.col(score,ties.method="first")]
tcell <- data.frame(Celltype=ct,UMAP1=x,UMAP2=y)
label_t <- data.frame(Celltype=names(cols),
  x=c(1.1,.45,4.25,.4,-3.2,-4.95,-.45,-.2,-6.6,4.95),
  y=c(-1.38,1.39,2.48,-2.56,-4.0,-.20,-.22,3.17,1.59,-.38))
write.csv(tcell,file.path(res,"23_tcell_visual_example.csv"),row.names=FALSE)
write.csv(label_t,file.path(res,"23_tcell_label_positions.csv"),row.names=FALSE)

## 2. 全细胞图：恢复被遗漏的图；比例、坐标均是新生成的合成示例 ====
cols_all <- c("Tumor cells"="#0d0887","Macrophages/Monocytes"="#33a02c","Neutrophils"="#88b817",
  "Stromal cells"="#0067a5","T cells"="#be0032","Vasculature"="#e68fac","Low count/Necrotic"="#d3d3d3",
  "Dendritic cells"="#fe00fe","NK cells"="#f48400","Plasma cells"="#00a9b5","B cells"="#f3c300")
# 每类轮廓和密度是图形参数，不能作为肿瘤组织/细胞比例估计。
regions <- list(
  "Tumor cells"=rbind(c(-3.3,1.5),c(-2.4,3.5),c(-.8,4.5),c(1.6,5.9),c(3.6,4.7),c(6.2,4.0),c(8.7,3.4),c(8.3,1.5),c(9.7,-.2),c(9.1,-1.9),c(7.0,-2.4),c(4.2,-4.5),c(1.3,-5.3),c(-.7,-4.5),c(-1.5,-2.9),c(-3.2,-2.3)),
  "Macrophages/Monocytes"=rbind(c(-6.5,-.4),c(-5.3,.7),c(-3.6,.8),c(-2.7,-.8),c(-1.2,-1.8),c(-2.5,-3.6),c(-4.7,-3.7),c(-6.2,-2.6)),
  "Neutrophils"=rbind(c(-6.1,1.0),c(-4.9,1.9),c(-3.2,2.7),c(-2.0,.8),c(-2.9,-1.5),c(-4.4,-2.5),c(-5.6,-1.2)),
  "Stromal cells"=rbind(c(-9.7,-1.0),c(-8.4,-1.2),c(-6.9,-2.6),c(-6.0,-3.8),c(-7.6,-4.6),c(-8.7,-3.8)),
  "T cells"=rbind(c(-10,2.9),c(-9.2,4.6),c(-7.4,5.1),c(-5.8,4.8),c(-6.3,2.5),c(-6.8,.7),c(-8.6,.8),c(-9.9,1.2)),
  "Vasculature"=rbind(c(-8.8,-2.2),c(-7.3,-3.1),c(-7.6,-4.9),c(-5.2,-5.6),c(-4.1,-6.5),c(-5.7,-7.3),c(-7.6,-6.9),c(-8.1,-5.0)),
  "Low count/Necrotic"=rbind(c(-4.9,1.3),c(-4.7,3.6),c(-3.6,3.7),c(-2.8,2.6),c(-1.3,2.8),c(-.4,1.3),c(-1.7,-.2),c(-3.5,.3)),
  "Dendritic cells"=rbind(c(-8.5,1.4),c(-6.7,2.5),c(-5.1,1.9),c(-4.4,.1),c(-6.1,-.5),c(-7.2,-.2)),
  "NK cells"=rbind(c(-9.6,2.7),c(-8.7,4.7),c(-7.1,5.7),c(-5.4,5.5),c(-5.5,4.1),c(-6.4,2.4),c(-7.3,1.8)),
  "Plasma cells"=rbind(c(-11.3,.6),c(-10.0,1.7),c(-8.7,.4),c(-7.3,-.1),c(-6.4,-1.5),c(-7.9,-2.6),c(-9.3,-1.2),c(-10.4,-.9)),
  "B cells"=rbind(c(-11.5,.6),c(-10.7,1.6),c(-9.7,2.0),c(-8.6,1.5),c(-8.2,.9),c(-9.5,-.2),c(-10.4,.1)))
n_all <- 4*c(65000,12000,6000,3500,12000,6500,6000,2800,800,3600,2500)
all_list <- list()
for(k in seq_along(regions)) {
  pp <- regions[[k]];tt <- seq(0,1,length.out=nrow(pp)+1)
  o <- list(x=spline(tt,c(pp[,1],pp[1,1]),method="periodic",n=250)$y,
            y=spline(tt,c(pp[,2],pp[1,2]),method="periodic",n=250)$y)
  n <- n_all[k];ncan <- n*8
  xx <- runif(ncan,min(pp[,1]),max(pp[,1]));yy <- runif(ncan,min(pp[,2]),max(pp[,2]))
  take <- which(sp::point.in.polygon(xx,yy,o$x,o$y)>0)[seq_len(n)]
  all_list[[k]] <- data.frame(Celltype=names(regions)[k],UMAP1=xx[take]+rnorm(n,sd=.06),UMAP2=yy[take]+rnorm(n,sd=.06))
}
all_df <- bind_rows(all_list)
label_all <- data.frame(Celltype=names(cols_all),x=c(4.0,-2.8,-2.55,-9.0,-9.0,-7.4,-2.5,-6.5,-6.6,-8.75,-10.8),
  y=c(.6,-1.65,.08,-2.75,2.9,-6.5,2.4,1.0,4.25,-.9,1.0))
write.csv(all_df,file.path(res,"23_all_cell_visual_example.csv"),row.names=FALSE)
write.csv(label_all,file.path(res,"23_all_cell_label_positions.csv"),row.names=FALSE)

## 3. 点层独立栅格化，标签/坐标仍为矢量；每个源图输出独立 PNG+PDF ====
for(view in c("all","tcell")) {
  dd <- if(view=="all") all_df else tcell
  lab <- if(view=="all") label_all else label_t
  pal <- if(view=="all") cols_all else cols
  xr <- if(view=="all") c(-12.5,10.5) else c(-8.8,7.1)
  yr <- if(view=="all") c(-8,6.5) else c(-5.25,4.85)
  width <- if(view=="all") 220 else 199; height <- if(view=="all") 145 else 142
  # 捕获仅含散点的绘图区。背景透明，300dpi，PDF文字不被一起栅格化。
  capture <- ragg::agg_capture(width=ceiling(width/25.4*300),height=ceiling(height/25.4*300),
    res=300,background="transparent")
  par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
  plot.new();plot.window(xlim=xr,ylim=yr)
  set.seed(260123); ii <- sample(seq_len(nrow(dd)))
  ptcol <- scales::alpha(pal[dd$Celltype[ii]],if(view=="all") .25 else .96)
  points(dd$UMAP1[ii],dd$UMAP2[ii],pch=16,cex=if(view=="all") .19 else .16,col=ptcol)
  raster <- capture();dev.off()
  p <- ggplot()+annotation_raster(raster,xmin=xr[1],xmax=xr[2],ymin=yr[1],ymax=yr[2],interpolate=FALSE)+
    geom_label(data=lab,aes(x,y,label=Celltype,fill=Celltype),size=if(view=="all") 6.0 else 5.0,label.padding=unit(.65,"mm"),
      label.r=unit(.7,"mm"),linewidth=.32,colour="black")+
    scale_fill_manual(values=pal)+
    scale_x_continuous(limits=xr,breaks=if(view=="all") c(-10,-5,0,5,10) else c(-4,0,4),expand=c(0,0))+
    scale_y_continuous(limits=yr,breaks=if(view=="all") c(-8,-4,0,4) else c(-5,-2.5,0,2.5),expand=c(0,0))+
    coord_fixed()+labs(x="umap_1",y="umap_2")+
    theme_classic(base_size=12,base_family="Arial")+
    theme(legend.position="none",axis.text=element_text(size=13.5,colour="black"),
      axis.title=element_text(size=17),axis.line=element_line(linewidth=.55),
      axis.ticks=element_line(linewidth=.5),plot.margin=margin(2,3,2,2))
  stem <- if(view=="all") "23_all_celltypes_umap" else "23_tcell_subtype_umap"
  save_plot(p,file.path(fig,stem),width,height)
}
cat("VISUAL_REPAIR: full-cell and T-cell plots;",nrow(tcell),"+",nrow(all_df),"synthetic points.\n")

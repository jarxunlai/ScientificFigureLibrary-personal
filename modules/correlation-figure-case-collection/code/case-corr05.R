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
# 第二轮：移除上缘锯齿、恢复模块三角边界、树标题和来源负相关色阶；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure/case-corr05")
dat <- file.path(module_root,"data/case-corr05")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 对数丰度的相关性与complete-linkage聚类 ====
d <- read.csv(file.path(dat,"abundance.csv"),row.names=1,check.names=FALSE)
g <- read.csv(file.path(dat,"groups.csv"))
r <- cor(log10(d+1e-6)); hc <- hclust(as.dist(1-r),method="complete"); ord <- hc$order; r <- r[ord,ord]; n <- nrow(r)
cells <- expand.grid(i=1:n,j=1:n); cells <- cells[cells$i<cells$j,]
cells$x <- (cells$i+cells$j)/2; cells$y <- -(cells$j-cells$i)/2
cells$value <- r[cbind(cells$i,cells$j)]
polys <- do.call(rbind,lapply(seq_len(nrow(cells)),function(k) data.frame(id=k,x=cells$x[k]+c(-.5,0,.5,0),y=cells$y[k]+c(0,.5,0,-.5),r=cells$value[k])))
tree <- ggdendro::dendro_data(as.dendrogram(hc))$segments
tree$y <- .8+tree$y/max(tree$y)*3; tree$yend <- .8+tree$yend/max(ggdendro::dendro_data(as.dendrogram(hc))$segments$y)*3
lab <- data.frame(k=1:n,name=colnames(r),group=g$group[match(colnames(r),g$variable)])
## 2. 用真实菱形多边形旋转矩阵，不旋转文字或贴参考图片 ====
p <- ggplot()+geom_polygon(data=polys,aes(x,y,group=id,fill=r),colour=NA)+
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=-.49,ymax=.60,fill="white",colour=NA)+
  geom_segment(data=tree,aes(x,y,xend=xend,yend=yend),linewidth=.27,colour="grey35")+
  geom_text(data=lab,aes(x=(k+1)/2-.65,y=-(k-1)/2-.45,label=name),angle=45,hjust=1,size=2.45)+
  geom_text(data=lab,aes(x=(k+n)/2+.65,y=-(n-k)/2-.45,label=name),angle=-45,hjust=0,size=2.45)+
  geom_point(data=lab,aes(x=(k+1)/2-.3,y=-(k-1)/2-.3,colour=group),size=1)+
  geom_point(data=lab,aes(x=(k+n)/2+.3,y=-(n-k)/2-.3,colour=group),size=1)+
  scale_fill_gradient2(low="#bbdeec",mid="white",high="#802520",limits=c(-.2,1),oob=scales::squish,name="Pearson r")+
  scale_colour_manual(values=c("#A69900","#8F1D00","#056258","#6A478F"),name="Variable group")+
  coord_fixed(xlim=c(-3,n+3),ylim=c(-n/2-5,5),clip="off")+
  theme_void(base_size=8)+theme(legend.position="right",plot.margin=margin(5,12,5,12))
modules<-cutree(hc,k=4)[ord]
for(k in unique(modules)) {ids<-which(modules==k);a<-min(ids)-.5;b<-max(ids)+.5
 p<-p+annotate("path",x=c(a,(a+b)/2,b),y=c(-.5,-(b-a)/2,-.5),linewidth=.4,colour="#333333")}
p<-p+annotate("text",x=n/2,y=4.4,label="Cluster Dendrogram",size=3.5)
save_plot(p,file.path(fig,"05_inverted_triangle_dendrogram"),220,155)
write.csv(data.frame(variable=colnames(r),position=1:n),file.path(root,"result/case-corr05","cluster_order.csv"),row.names=FALSE)

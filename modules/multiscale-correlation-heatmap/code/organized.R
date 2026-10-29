# 一维二维多维相关性热图 / Multiscale correlation heatmaps
# 科学问题：单组变量、两类特征及群落与环境的关联分别如何展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "multiscale-correlation-heatmap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：补齐右侧分组注释、显著性留白和完整色标，恢复聚类树与Mantel图例可读性；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 一维：按变量组着色的显著相关三角矩阵 ====
d <- read.csv(file.path(dat,"environment_variables.csv"),check.names=FALSE)
groups <- read.csv(file.path(dat,"environment_groups.csv"))
m <- scale(as.matrix(d[,-1]));set.seed(270102)
# 降低原示例过强的共同因子，避免显著性矩阵整块发黑。
m <- .82*m+matrix(rnorm(length(m),sd=.60),nrow(m),ncol(m));r <- cor(m)
write.csv(m,file.path(entry,"result/02a_visual_observations.csv"))
# 文章把相关矩阵再次送入cor.mtest；此处直接对样本观测列检验，避免二次相关。
p <- matrix(0,ncol(m),ncol(m),dimnames=dimnames(r))
for(i in 1:ncol(m)) for(j in 1:ncol(m)) p[i,j] <- cor.test(m[,i],m[,j])$p.value
group_cols <- c(Microbial="#66009b",Site="#00bccc",`Land use`="#ff9933",Watershed="#00089d")
save_base({
 par(mar=c(0,0,0,0),ps=11,xaxs="i",yaxs="i")
 n<-ncol(r);plot.new();plot.window(xlim=c(-16,n+20),ylim=c(-8,n+18),asp=1)
 pal<-circlize::colorRamp2(c(-1,0,1),c("#b2182b","white","#1a1a1a"))
 for(i in 1:n) for(j in i:n){
   rect(j-.5,n-i+.5,j+.5,n-i+1.5,col="white",border="#cccccc",lwd=.3)
   if(p[i,j]<.01){h<-.46*sqrt(abs(r[i,j]));rect(j-h,n-i+1-h,j+h,n-i+1+h,col=pal(r[i,j]),border=NA)}
 }
 for(i in 1:n){
   text(i-.65,n-i+1,colnames(r)[i],adj=1,cex=.59,col=group_cols[groups$group[i]])
   text(i,n+.8,colnames(r)[i],srt=90,adj=0,cex=.60,col=group_cols[groups$group[i]])
   rect(n+1,n-i+.5,n+2.1,n-i+1.5,col=group_cols[groups$group[i]],border=NA)
   text(n+2.7,n-i+1,colnames(r)[i],adj=0,cex=.59)
 }
 z<-seq(-1,1,length.out=101)
 rect(4+head(z+1,-1)*n*.25,-4,4+tail(z+1,-1)*n*.25,-2.8,col=pal(head(z,-1)),border=NA)
 text(4+(c(-1,-.5,0,.5,1)+1)*n*.25,-5.5,c(-1,-.5,0,.5,1),cex=.65)
},file.path(fig,"02a_grouped_environment_correlation"),200,185)
## 2. 二维：NMF因子与细胞类型，使用对称颜色范围 ====
cell <- read.csv(file.path(dat,"cell_factor_correlation.csv"),row.names=1,check.names=FALSE)
annotation <- read.csv(file.path(dat,"cell_groups.csv"),row.names=1)
pal <- colorRampPalette(c("#123d35","#72b5a1","white","#d296c1","#753a75"))(101)
lim <- max(abs(as.matrix(cell)))
hm <- pheatmap::pheatmap(as.matrix(cell),color=pal,breaks=seq(-lim,lim,length.out=102),annotation_row=annotation,
  annotation_colors=list(lineage=c(Epithelial="#fbb03c",Endothelial="#a0e0b8",Immune="#0070bc",Mesenchymal="#b36692",Other="#faee90")),
  annotation_names_row=FALSE,treeheight_row=18,treeheight_col=18,border_color="grey85",fontsize=9,angle_col=90,silent=TRUE)
save_base({grid.newpage(); grid.draw(hm$gtable)},file.path(fig,"02b_cell_factor_clustered_heatmap"),190,160)
## 3. 多维：控制合成地理坐标的partial Mantel；非真实群落上游分析 ====
env <- read.csv(file.path(dat,"mantel_environment.csv"),row.names=1,check.names=FALSE)
spec <- read.csv(file.path(dat,"mantel_communities.csv"),row.names=1,check.names=FALSE)
geo <- read.csv(file.path(dat,"mantel_geography.csv"),row.names=1)
set.seed(260102)
mt <- linkET::mantel_test(spec,env,env_ctrl=geo,mantel_fun="mantel.partial",spec_select=list(host=1:6,mOTUs=7:12,virus=13:18,vOTUs=19:24),permutations=999)
mt$rd <- cut(mt$r,breaks=c(-Inf,.1,.2,Inf),right=FALSE,labels=c("< 0.1","[0.1, 0.2)",">= 0.2"))
# 999次置换的最小P为0.001；默认cut右闭区间必须明确写成<=，不能标成<。
mt$pd <- cut(mt$p,breaks=c(-Inf,.001,.05,Inf),labels=c("<= 0.001","(0.001, 0.05]","> 0.05"))
write.csv(mt,file.path(root,"result","mantel_results.csv"),row.names=FALSE)
p3 <- linkET::qcorrplot(linkET::correlate(env),type="upper",diag=FALSE)+
  linkET::geom_square()+linkET::geom_mark(only_mark=TRUE,size=1.7)+
  linkET::geom_couple(aes(colour=pd,size=rd),data=mt,curvature=.12,node.colour="grey30",node.fill="white",node.size=2)+
  scale_fill_gradient2(low="#ff8040",mid="white",high="#5bc2cd",limits=c(-1,1),name="Pearson r")+
  scale_size_manual(values=c(.35,.8,1.4),name="Mantel r")+
  scale_colour_manual(values=c("#87ceeb","#9acd32","#cfcfcf"),name="Partial Mantel P")+
  theme(legend.title=element_text(size=9),legend.text=element_text(size=8),axis.text=element_text(size=8))
save_plot(p3,file.path(fig,"02c_mantel_environment_network"),215,160)

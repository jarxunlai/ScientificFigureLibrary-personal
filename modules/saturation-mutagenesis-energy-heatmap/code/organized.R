# 饱和突变能量热图与均值折线 / Saturation mutagenesis energy heatmap with site means
# 科学问题：哪些残基位点的替换总体更敏感，哪些具体替换具有较高预测能量代价？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "saturation-mutagenesis-energy-heatmap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复均值正负点和零线的清晰层级、氨基酸/位点标签及双轴单位；输入/派生值均为合成示例，不冒充来源结果。
# 饱和突变能量热图：20 AA × 64 位点，均值折线来自同一输入矩阵。
# 网页 heat_fine 未给定义；本稿显式双线性插值用于背景显示，统计量始终从原始离散矩阵算。
# 无分子模拟、对接或真实突变实验；独立模拟能量输入用于绘图模板。
library(ggplot2)
library(grid)


source(file.path(module_root,"code/export.R"))

## 1. 输入、按位点均值与可逆双轴映射 ====
d<-read.csv(file.path(module_root,"data/mutation_energy.csv"),check.names=FALSE)
m<-as.matrix(d[,-1]);rownames(m)<-d$MutAA
nx<-ncol(m);ny<-nrow(m)
means<-data.frame(x=1:nx,Site=colnames(m),MeanDDG=colMeans(m))
ddg_low<-floor(min(c(means$MeanDDG,0))*10)/10
ddg_high<-ceiling(max(c(means$MeanDDG,0))*10)/10
means$y<-1+(means$MeanDDG-ddg_low)/(ddg_high-ddg_low)*(ny-1)
y_zero<-1+(0-ddg_low)/(ddg_high-ddg_low)*(ny-1)
write.csv(means,file.path(draft,"result/11_site_means.csv"),row.names=FALSE)
write.csv(data.frame(AA=d$MutAA,m,check.names=FALSE),file.path(draft,"result/11_discrete_plot_input.csv"),row.names=FALSE)

## 2. 显式插值：先沿位点再沿氨基酸轴；不改变均值或新增统计 ====
x_fine<-seq(1,nx,length.out=nx*6)
y_fine<-seq(1,ny,length.out=ny*6)
step_x<-matrix(0,ny,length(x_fine))
for(i in 1:ny) step_x[i,]<-approx(1:nx,m[i,],xout=x_fine)$y
fine<-matrix(0,length(y_fine),length(x_fine))
for(j in seq_along(x_fine)) fine[,j]<-approx(1:ny,step_x[,j],xout=y_fine)$y
heat_fine<-expand.grid(y=y_fine,x=x_fine)
heat_fine$DDG<-as.vector(fine)
cols24<-c("#008100","#0d930f","#2fa721","#40c52e","#59dc41","#6aed51","#82fe62","#99fe72","#b4fb83","#ccfe90","#dbfca1","#f6fcb2","#fef5b0","#fddda4","#ffcd91","#feb483","#ff9a70","#ff8160","#fd6852","#fe5841","#fe4030","#ff2720","#fd0f0f","#f80103")

## 3. 热图、零线、灰色均值折线与粉/蓝点 ====
p<-ggplot()+
  geom_raster(data=heat_fine,aes(x,y,fill=DDG),interpolate=TRUE)+
  scale_fill_gradientn(colours=cols24,limits=c(-1.5,3),breaks=seq(-1.5,3,.5),
                       name="Color scales of saturated mutation energy")+
  geom_hline(yintercept=y_zero,linetype="dashed",linewidth=.25,colour="black")+
  geom_line(data=means,aes(x,y,group=1),linetype="dashed",linewidth=.25,colour="#555555")+
  geom_point(data=means,aes(x,y,colour=MeanDDG>0),size=1.2)+
  scale_colour_manual(values=c("TRUE"="#f80660","FALSE"="#2c43ec"),guide="none")+
  scale_x_continuous(breaks=1:nx,labels=colnames(m),expand=c(0,0))+
  scale_y_continuous(breaks=1:ny,labels=rownames(m),expand=c(0,0),
      sec.axis=sec_axis(transform=~(. - 1)/(ny-1)*(ddg_high-ddg_low)+ddg_low,
                        name="Average mutation energy\n(binding, kcal/mol)",breaks=pretty(c(ddg_low,ddg_high),n=6)))+
  coord_cartesian(xlim=c(.5,nx+.5),ylim=c(.5,ny+.5),expand=FALSE)+
  labs(x=NULL,y="Mutation sites")+
  theme_classic(base_family="Arial",base_size=8)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8.2,colour="black"),
        axis.text.y=element_text(size=7,colour="black"),axis.title=element_text(size=8),
        axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),axis.ticks.length=unit(.8,"mm"),
        legend.position="top",legend.title=element_text(size=8,face="bold"),legend.text=element_text(size=7),
        plot.margin=margin(2,2,2,2,"mm"))+
  guides(fill=guide_colourbar(direction="horizontal",title.position="top",barwidth=unit(145,"mm"),
                             barheight=unit(3,"mm"),ticks=TRUE))
save_gg(p,file.path(draft,"figure/11_mutation_energy"),180,99)

# 四种散点图局部放大布局 / Four local-magnification scatter layouts
# 科学问题：如何同时保留全局分布并显示高密度区域中的关键细节？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "local-magnification-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：补回古代PCA人群图例，修复外置放大窗画布留白和四种布局的标签可读性；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 局部放大 - 火山图/PCA图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1025)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df<-read.csv(file.path(module_root,"data/dependency.csv"))
df$sig <- df$diff<=-.3 & df$score>=-log10(.05)
## 3. 案例一：右侧火山图放大框；源窗口与引线显式绘制 ====
pm <- ggplot(df,aes(diff,score)) + geom_point(data=subset(df,!sig),colour="grey60",size=.3,alpha=.3) + geom_point(data=subset(df,sig),colour="#a6447c",size=1,alpha=.75) +
  geom_hline(yintercept=-log10(.05),linetype=3,linewidth=.3) + geom_vline(xintercept=-.3,linetype=3,linewidth=.3) +
  geom_text_repel(data=subset(df,label_main==1),aes(label=gene),size=2.3,colour="#a6447c",seed=65,max.overlaps=Inf) +
  annotate("rect",xmin=-.56,xmax=-.3,ymin=1.3,ymax=5.5,fill=NA,colour="#ed5c79",linetype=3,linewidth=.35) + coord_cartesian(xlim=c(-1,1),ylim=c(0,15)) + labs(x="Delta dependency (glial - other NextGen)",y="-log10 P")
zoom_data <- subset(df,diff>=-.56 & diff<=-.3 & score>=1.3 & score<=5.5)
# 显式标注前30个示例标记；所有窗口内的点仍显示，不对密集点逐一强加标签。
pz <- ggplot(zoom_data,aes(diff,score)) + geom_point(colour="#a6447c",size=1.1,alpha=.6) + geom_text_repel(data=subset(zoom_data,zoom_label==1),aes(label=gene),size=2.5,seed=65,max.overlaps=Inf,min.segment.length=0) +
  coord_cartesian(xlim=c(-.56,-.3),ylim=c(1.3,5.5),expand=FALSE) + labs(x=NULL,y=NULL) + theme(panel.border=element_rect(fill=NA,colour="#ed5c79",linewidth=.4),axis.text=element_text(size=5))
pc <- pm + coord_cartesian(xlim=c(-1,1),ylim=c(0,15),clip="off",expand=FALSE) + theme(plot.margin=margin(8,75,8,5,unit="mm")) +
  annotation_custom(ggplotGrob(pz),xmin=1.2,xmax=2.75,ymin=1,ymax=14) +
  annotate("segment",x=-.3,xend=1.2,y=5.5,yend=14,colour="#ed5c79",linetype=3,linewidth=.3) +
  annotate("segment",x=-.3,xend=1.2,y=1.3,yend=1,colour="#ed5c79",linetype=3,linewidth=.3)
save_plot(pm,file.path(fig_dir,"25a_dependency_main"),100,100);save_plot(pz,file.path(fig_dir,"25a_dependency_zoom"),70,100);save_plot(pc,file.path(fig_dir,"25a_dependency_magnified"),170,105)
## 4. 案例二：PCA中心样本左侧放大（输入为预计算scores） ====
dd<-read.csv(file.path(module_root,"data/pca_samples.csv"));pal<-setNames(scales::hue_pal()(10),unique(dd$Group))
pm<-ggplot(dd,aes(PC1,PC2,colour=Group)) + geom_point(size=1.3,alpha=.8) + scale_colour_manual(values=pal) +
  annotate("rect",xmin=-17,xmax=20.5,ymin=-7,ymax=20,fill=NA,colour="grey30",linetype=3) + labs(x="PC1",y="PC2") + theme(legend.position="bottom",legend.text=element_text(size=5),legend.key.size=unit(2,"mm"))
pz<-ggplot(subset(dd,PC1>=-17 & PC1<=20.5 & PC2>=-7 & PC2<=20),aes(PC1,PC2,colour=Group)) + geom_point(size=2.5) + scale_colour_manual(values=pal) + theme_void() + theme(legend.position="none",panel.border=element_rect(fill=NA,colour="grey35",linewidth=.3))
pc<-pm + coord_cartesian(xlim=c(-50,80),ylim=c(-45,45),clip="off",expand=FALSE) + theme(plot.margin=margin(5,5,10,72,unit="mm")) +
  annotation_custom(ggplotGrob(pz),xmin=-150,xmax=-60,ymin=-25,ymax=35) +
  annotate("segment",x=-17,xend=-60,y=20,yend=35,colour="grey50",linetype=3,linewidth=.3) + annotate("segment",x=-17,xend=-60,y=-7,yend=-25,colour="grey50",linetype=3,linewidth=.3)
save_plot(pc,file.path(fig_dir,"25b_pca_left_magnifier"),180,100)
## 5. 案例三：古代DNA PCA中突出LBA样本与放大窗口 ====
anc<-read.csv(file.path(module_root,"data/ancient_pca.csv"));pal2<-setNames(scales::hue_pal()(15),unique(anc$Population))
pm<-ggplot(anc,aes(PC1,PC2)) + geom_point(aes(colour=Population,shape=Population),size=1.2,alpha=.35) +
  geom_point(data=subset(anc,New),aes(fill=Population),shape=21,size=2,colour="black",stroke=.4) + scale_colour_manual(values=pal2) + scale_fill_manual(values=pal2) + scale_shape_manual(values=rep(c(16,17,15,18,1),3)) +
  annotate("rect",xmin=-.035,xmax=-.014,ymin=-.027,ymax=.005,fill=NA,colour="grey35",linewidth=.3) + labs(x="PC1",y="PC2") + theme(legend.position="bottom",legend.text=element_text(size=5.8),legend.key.size=unit(2,"mm"))
zz<-subset(anc,PC1>=-.035 & PC1<=-.014 & PC2>=-.027 & PC2<=.005)
pz<-ggplot(zz,aes(PC1,PC2)) + geom_point(aes(colour=Population),size=1,alpha=.25) + geom_point(data=subset(zz,New),aes(fill=Population),shape=21,size=2,colour="black",stroke=.5) +
  scale_colour_manual(values=pal2) + scale_fill_manual(values=pal2) + geom_text_repel(data=head(subset(zz,New),3),aes(label=ID),size=2.3,seed=66) + theme_void() + theme(legend.position="none",panel.border=element_rect(fill=NA,colour="grey35",linewidth=.3))
pc<-pm + coord_cartesian(xlim=c(-.14,.08),ylim=c(-.09,.09),clip="off",expand=FALSE) + theme(plot.margin=margin(5,75,22,5,unit="mm")) +
  annotation_custom(ggplotGrob(pz),xmin=.1,xmax=.25,ymin=-.075,ymax=.075) +
  annotate("segment",x=-.014,xend=.1,y=.005,yend=.075,colour="grey50",linewidth=.3) + annotate("segment",x=-.014,xend=.1,y=-.027,yend=-.075,colour="grey50",linewidth=.3)
save_plot(pc,file.path(fig_dir,"25c_ancient_pca_right_magnifier"),180,110)
## 6. 案例四：极端值总览在上，密集区放大在下 ====
dd<-read.csv(file.path(module_root,"data/neuropeptides.csv"))
pm<-ggplot(dd,aes(Neuropeptides,Receptors)) + geom_point(colour="#5d3a9b",size=1) + geom_text_repel(data=subset(dd,extreme),aes(label=ID),size=2.3,seed=67) +
  annotate("rect",xmin=0,xmax=25,ymin=0,ymax=50,fill=NA,colour="grey40",linetype=3,linewidth=.3) + coord_cartesian(xlim=c(0,100),ylim=c(0,80)) + labs(x=NULL,y="Cognate receptors\n(% NE cells)") + theme(axis.title.y=element_text(size=6))
pz<-ggplot(subset(dd,!extreme),aes(Neuropeptides,Receptors)) + geom_point(colour="#5d3a9b",size=1.8) + geom_text_repel(data=subset(dd,!extreme & labeled),aes(label=ID),size=2.2,fontface="italic",seed=67,max.overlaps=Inf,min.segment.length=0) +
  coord_cartesian(xlim=c(0,25),ylim=c(0,50)) + labs(x="Neuropeptides (% GABA cells)",y="Cognate receptors (% NE cells)") + theme(axis.title=element_text(size=7))
pc<-cowplot::ggdraw() + cowplot::draw_plot(pm,.25,.60,.53,.38) + cowplot::draw_plot(pz,0,0,1,.58) + cowplot::draw_line(x=c(.35,.10),y=c(.68,.58),colour="grey50",linetype=3,linewidth=.3) + cowplot::draw_line(x=c(.45,.90),y=c(.68,.58),colour="grey50",linetype=3,linewidth=.3)
save_plot(pm,file.path(fig_dir,"25d_neuropeptide_overview"),85,70);save_plot(pz,file.path(fig_dir,"25d_neuropeptide_zoom"),130,90);save_plot(pc,file.path(fig_dir,"25d_overview_and_zoom"),140,150)

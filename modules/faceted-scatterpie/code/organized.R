# 医院与质粒类型分面散点饼图 / Faceted hospital-plasmid scatter pies
# 科学问题：耐药基因、医院、复制子质粒和菌种组成如何联合分布？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "faceted-scatterpie")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复菌株数气泡半径及对应图例，放大医院标签，保留三条独立耐药基因分面；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 分面散点饼图全网最详细教程（附代码，文末有福利）
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1008)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

library(scatterpie)
df <- read.csv(file.path(module_root,"data/hospital_plasmids.csv"),check.names=FALSE)
species <- c("Klebsiella spp","Escherichia coli","Enterobacter spp","Citrobacter freundii","others")
cols <- c("green3","darkorange1","deeppink","blue2","grey40")
panels <- list()
## 3. 各耐药基因独立y轴：气泡半径随菌株数平方根变化 ====
for (gene in c("blaIMP-1","blaIMP-6","blaNDM-5")) {
  dd <- subset(df,Gene==gene);dd$Radius<-dd$Radius*1.55
  key<-data.frame(x0=47,r=1.55*(.10+.24*sqrt(c(5,15,30,45)/max(df$Total))),label=c(5,15,30,45))
  key$y0<-.55+key$r
  p <- ggplot() + ggforce::geom_circle(data=key,aes(x0=x0,y0=y0,r=r),colour="grey40",fill=NA,linewidth=.15)+
    geom_text(data=key,aes(x=47.8,y=.55+2*r,label=label),hjust=0,size=1.7)+geom_scatterpie(data=dd,aes(x=Hospital_id,y=Plasmid_id,r=Radius),cols=species,colour="grey30",linewidth=.12,alpha=.8) +
    scale_fill_manual(values=setNames(cols,species)) + scale_x_continuous(breaks=1:45,labels=sprintf("Hospital %02d",1:45),limits=c(.3,49.5)) +
    scale_y_continuous(breaks=sort(unique(dd$Plasmid_id)),labels=unique(dd$Plasmid[order(dd$Plasmid_id)])) +
    labs(y=gene,x=NULL) + coord_fixed() + theme_bw(base_size=7) + theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour="grey80",linewidth=.2),axis.text.x=element_text(angle=90,hjust=1,size=5.8),legend.text=element_text(size=6))
  panels[[gene]] <- p+theme(axis.text.x=if(gene=="blaNDM-5")element_text(angle=90,hjust=1,size=5.8)else element_blank(),axis.ticks.x=if(gene=="blaNDM-5")element_line()else element_blank())
  save_plot(p,file.path(fig_dir,paste0("08_",gene)),180,60)
}
radius_values <- c(5,15,30,45)
size_legend <- data.frame(x=1:4,y=1,r=1.55*(.10+.24*sqrt(radius_values/max(df$Total))),label=radius_values)
size_plot <- ggplot(size_legend,aes(x,y)) + geom_point(aes(size=r*5),shape=1,colour="grey40",stroke=.25) + scale_size_identity() +
  geom_text(aes(y=.6,label=label),size=2) + coord_cartesian(xlim=c(.5,4.5),ylim=c(.4,1.5),expand=FALSE) + labs(title="Number of strains") + theme_void() + theme(plot.title=element_text(size=6,hjust=.5))
save_plot(size_plot,file.path(fig_dir,"08_strain_size_legend"),75,30)
species_key<-data.frame(Species=species,y=5:1)
species_plot<-ggplot(species_key,aes(0,y)) + geom_point(aes(fill=Species),shape=22,size=2.5,colour="grey50",stroke=.2) +
  geom_text(aes(x=.18,label=Species),hjust=0,size=1.9) + scale_fill_manual(values=setNames(cols,species),guide="none") + xlim(-.1,1.5) + theme_void()
pie_composite <- wrap_plots(panels,ncol=1,heights=c(9,4,4)) & theme(legend.position="none")
legend_column <- species_plot + labs(caption="Circle key: number of strains")+theme(plot.caption=element_text(size=5))
save_plot(wrap_plots(pie_composite,legend_column,ncol=2,widths=c(5,1.25)),file.path(fig_dir,"08_faceted_scatterpie"),205,85)

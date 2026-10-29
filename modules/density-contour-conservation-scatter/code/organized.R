# 填充密度等高线与跨物种保守性 / Filled density contours and cross-species conservation
# 科学问题：二维变量的密集区域在哪里，不同保守性类别偏离对角线多少？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "density-contour-conservation-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：固定四类等高线顺序并恢复密度层数、保守性参考线和标记标签层级；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 等高线图也可以很炫酷！
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1013)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/periodicity.csv"))
panels <- list()
## 3. 四类区域的二维核密度等高线 ====
for (g in c("All","PLD","Aromatic","TFs")) {
  dd <- subset(df,Region==g)
  col <- c(All="#1466a1",PLD="#593b94",Aromatic="#d7191c",TFs="#238b45")[g]
  p <- ggplot(dd,aes(logP,length)) + stat_density_2d_filled(contour_var="ndensity",bins=12,show.legend=FALSE) +
    scale_fill_manual(values=colorRampPalette(c("white",col))(12)) +
    annotate("text",x=4.8,y=285,label=paste0("(",nrow(dd),")"),hjust=1,size=2.3) +
    coord_cartesian(xlim=c(2,5),ylim=c(100,300),expand=FALSE) + labs(title=g,x="-log10(P value)",y="Length (AA)") + theme(panel.border=element_rect(fill=NA,colour="grey55",linewidth=.3))
  panels[[g]] <- p; save_plot(p,file.path(fig_dir,paste0("13_",g)),55,70)
}
save_plot(wrap_plots(panels,nrow=1),file.path(fig_dir,"13a_periodicity_contours"),180,72)
## 4. 跨物种保守性：灰色轮廓、标记基因及对角参考 ====
dd <- read.csv(file.path(module_root,"data/conservation.csv"))
p <- ggplot(dd,aes(human,lemur)) + stat_density_2d(aes(colour=after_stat(level)),bins=16,linewidth=.28) +
  scale_colour_gradient(low="#bab9b0",high="black",name="Contour density") + ggnewscale::new_scale_colour() +
  geom_abline(slope=1,intercept=c(-.4,0,.4),linetype=c(3,2,3),colour=c("#d66a61","grey40","#6c7bb8"),linewidth=.3) +
  geom_point(data=subset(dd,category!="background"),aes(colour=category),size=.75) +
  geom_text_repel(data=subset(dd,label==1),aes(label=gene,colour=category),size=2.4,seed=53,max.overlaps=Inf,min.segment.length=0) +
  scale_colour_manual(values=c("Highly conserved"="#20b521","HL conserved"="#d66a61","HM conserved"="#6c7bb8")) + coord_fixed(xlim=c(-.3,1),ylim=c(-.3,1)) +
  labs(x="rg(HL)",y="rg(HM)") + theme(legend.position="right",legend.text=element_text(size=6))
save_plot(p,file.path(fig_dir,"13b_conservation_contours"),170,130)

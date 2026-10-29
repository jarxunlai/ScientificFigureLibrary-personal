# 配对散点与对角线旋转直方图 / Paired scatter with diagonal rotated histograms
# 科学问题：处理前后配对读数的差异在条件和动物之间如何变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "paired-scatter-rotated-histograms")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复六面板比例，放大对角旋转直方图及其中位数标记，保留空实心配对点；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 对角线散点图+直方图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1018)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/paired_neural.csv"))
cols <- c(Vestibular="#327abd",Visual="#f71d53",Combined="#3b9e64")
panels <- list()
## 3. 六个配对散点面板；对角线外放置旋转差值直方图 ====
for (animal in c("Monkey D","Monkey F")) for (condition in names(cols)) {
  dd <- subset(df,Animal==animal & Condition==condition)
  rr <- wilcox.test(dd$Treat,dd$Control,paired=TRUE,exact=FALSE)
  med <- median(dd$Treat-dd$Control);col <- cols[condition]
  pm <- ggplot(dd,aes(Control,Treat)) + geom_abline(slope=1,intercept=0,linetype=2,linewidth=.3) + geom_hline(yintercept=0,linetype=3,linewidth=.25) + geom_vline(xintercept=0,linetype=3,linewidth=.25) +
    geom_point(aes(fill=Significant),shape=if(animal=="Monkey D")21 else 24,colour=col,size=1.8,stroke=.4) + scale_fill_manual(values=c("no"="white","yes"=col)) +
    annotate("text",x=-.95,y=.95,label=sprintf("P = %.2g",rr$p.value),hjust=0,size=2.3) + coord_fixed(xlim=c(-1,1),ylim=c(-1,1)) +
    labs(x="PSE without microstim. (deg)",y="PSE with microstim. (deg)",title=paste(animal,condition,sep=" / ")) + theme(legend.position="none",axis.title=element_text(size=7),plot.title=element_text(size=7),plot.margin=margin(2,2,2,2))
  ph <- ggplot(dd,aes(Treat-Control,fill=Significant)) + geom_histogram(bins=10,position="stack",colour=col,linewidth=.25) + scale_fill_manual(values=c("no"="white","yes"=col)) +
    geom_vline(xintercept=med,linewidth=.4) + labs(x=NULL,y=NULL) + ggtitle(sprintf("median %.2f",med)) + theme_classic(base_size=7) + theme(legend.position="none",plot.title=element_text(size=7,hjust=.5),plot.margin=margin(0,0,0,0))
  rotated <- grobTree(ggplotGrob(ph),vp=viewport(angle=-45,width=.9,height=.85))
  panel <- cowplot::ggdraw() + cowplot::draw_plot(pm,x=0,y=0,width=.84,height=.84) + cowplot::draw_grob(rotated,x=.60,y=.64,width=.35,height=.35)
  key <- paste(gsub(" ","_",animal),condition,sep="_");panels[[key]]<-panel
  save_plot(panel,file.path(fig_dir,paste0("18_",key)),90,90)
}
save_plot(wrap_plots(panels,ncol=3),file.path(fig_dir,"18_neural_paired_six_panel"),210,140)

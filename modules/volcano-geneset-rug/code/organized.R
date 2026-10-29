# 火山图与基因集rug条 / Volcano plot with gene-set rug tracks
# 科学问题：差异表达分布与特征基因集成员在FC轴上如何对应？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "volcano-geneset-rug")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复KRAS目标标签和双侧Top区，压缩底部rug高度并保留共同横坐标；输入/派生值均为合成示例，不冒充来源结果。
# 如何把普通的火山图，做成顶刊里的样子？
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1009)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/volcano.csv"))
ix<-which(df$label==1);df$gene[ix]<-rep(c("KRAS","FOSL1","MYC"),length.out=length(ix))
rugs <- read.csv(file.path(module_root,"data/signature_rug.csv"))
cols <- c("KRAS signature"="#dfc956","Other"="#606060")
## 3. 火山图、双侧top区和两条特征基因rug ====
p <- ggplot(df,aes(logFC,score)) +
  annotate("rect",xmin=-3,xmax=-.7,ymin=1.3,ymax=6.5,fill="#caddec",alpha=.7,linetype=3,colour="grey60") +
  annotate("rect",xmin=.7,xmax=3,ymin=1.3,ymax=6.5,fill="#d8e6ac",alpha=.7,linetype=3,colour="grey60") +
  geom_point(aes(fill=signature),shape=21,colour="grey25",stroke=.15,size=.9,alpha=.65) +
  scale_fill_manual(values=cols) + geom_hline(yintercept=1.3,linetype=3,linewidth=.3,colour="grey60") +
  geom_label_repel(data=subset(df,label==1),aes(label=gene),size=2.6,label.padding=unit(.5,"mm"),min.segment.length=0,seed=49) +
  annotate("text",x=c(-1.8,1.8),y=6.8,label=c("PDAC KRAS UP\nTop 200","PDAC KRAS DOWN\nTop 200"),colour=c("#7ea9c6","#a9bc60"),size=2.8) +
  labs(x=NULL,y="Significance (-log10 adj. P)") + coord_cartesian(xlim=c(-3.4,3.4),ylim=c(0,7)) + theme(legend.position="none",axis.text.x=element_blank(),axis.ticks.x=element_blank())
pr <- ggplot(rugs,aes(logFC,Signature,colour=Signature)) + geom_point(shape="|",size=2.5,stroke=.3) +
  scale_colour_manual(values=c("KRAS_SIGNALING_UP"="#dfc956","KRAS_SIGNALING_DN"="#999999")) + labs(x="KRAS vs NS siRNA (log2FC)",y=NULL) +
  coord_cartesian(xlim=c(-3.4,3.4)) + theme_classic(base_size=7) + theme(legend.position="none",axis.text.y=element_text(size=5),axis.line.y=element_blank(),axis.ticks.y=element_blank())
save_plot(p,file.path(fig_dir,"09a_volcano"),120,100)
save_plot(pr,file.path(fig_dir,"09b_signature_rug"),120,25)
save_plot(p/pr+plot_layout(heights=c(6.5,.65)),file.path(fig_dir,"09_volcano_signature"),125,125)

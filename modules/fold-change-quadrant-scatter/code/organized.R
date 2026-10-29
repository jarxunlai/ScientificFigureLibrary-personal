# 转录翻译FC象限与联合显著性 / Transcription-translation fold-change quadrants
# 科学问题：转录与翻译的变化哪些一致、抵消或相反？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "fold-change-quadrant-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：修正六种变化模式计数框的上下分区，避免背景计数混入上调区域；保留象限语义；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | log2FC-FC 散点+象限图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1010)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/rna_fc.csv"))
cols <- c("Translation up"="#eb9778","Translation down"="#bd4d3c","Offsetting (RNA up)"="#8db7d8","Offsetting (RNA down)"="#4885b0","RNA abundance up"="#a4d18e","RNA abundance down"="#528e55","background"="grey70")
## 3. 转录/翻译/抵消分型：原始示例类别不是由图形反推 ====
p1 <- ggplot(df,aes(total,poly)) + geom_point(aes(colour=mode),size=.55,alpha=.6) +
  geom_hline(yintercept=c(-.5,.5),linewidth=.35) + geom_vline(xintercept=c(-.5,.5),linewidth=.35) + geom_abline(slope=1,intercept=c(-.5,.5),linewidth=.35) +
  geom_abline(slope=1,intercept=0,linetype=3,linewidth=.25) + geom_hline(yintercept=0,linetype=3,linewidth=.25) + geom_vline(xintercept=0,linetype=3,linewidth=.25) +
  geom_text_repel(data=subset(df,label==1),aes(label=gene),size=2,seed=50,max.overlaps=Inf,min.segment.length=0) + scale_colour_manual(values=cols) +
  coord_fixed(xlim=c(-5,5),ylim=c(-5,5)) + labs(title="shEif2b5 versus shCon",x="Total RNA (log2FC)",y="Polysome-associated RNA (log2FC)") + theme(legend.text=element_text(size=6))
counts <- table(factor(df$mode,levels=names(cols)))
cnt <- data.frame(mode=names(counts)[c(1,3,5,2,4,6)],n=as.integer(counts)[c(1,3,5,2,4,6)],x=-4.65,y=c(4.6,4.1,3.6,-3.6,-4.1,-4.6))
p1 <- p1 + geom_point(data=cnt,aes(x,y,fill=mode),inherit.aes=FALSE,shape=22,size=2,colour="black") + scale_fill_manual(values=cols,guide="none") +
  geom_text(data=cnt,aes(x=x+.35,y=y,label=n),inherit.aes=FALSE,hjust=0,size=2)
save_plot(p1,file.path(fig_dir,"10a_translation_quadrants"),160,132)
## 4. 独立FDR组合编码与一致/相反象限 ====
df$sig <- case_when(df$FDR_total<.05 & df$FDR_poly<.05~"both",df$FDR_total<.05~"totalRNA",df$FDR_poly<.05~"polyRNA",TRUE~"not significant")
p2 <- ggplot(df,aes(total,poly,colour=sig)) + geom_point(size=.6,alpha=.65) + geom_vline(xintercept=c(-.5,.5),linetype=3,linewidth=.4) + geom_hline(yintercept=c(-.5,.5),linetype=3,linewidth=.4) +
  scale_colour_manual(values=c(both="#a96cc6",polyRNA="#54b6bb",totalRNA="#cf97cb","not significant"="grey80")) +
  annotate("text",x=c(-3.6,3.5,-3.6,3.5),y=c(4.6,4.6,-4.6,-4.6),label=c("Discordant","Concordant","Concordant","Discordant"),size=2.5) +
  coord_fixed(xlim=c(-5,5),ylim=c(-5,5)) + labs(x="log2(Fold Change) totalRNA",y="log2(Fold Change) polyRNA") + theme(legend.text=element_text(size=6))
save_plot(p2,file.path(fig_dir,"10b_significance_quadrants"),150,128)
write.csv(df,file.path(root,"result/classified_points.csv"),row.names=FALSE)

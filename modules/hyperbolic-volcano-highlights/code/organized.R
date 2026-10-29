# 双曲线火山图与目标基因高亮 / Hyperbolic volcano plots and highlighted targets
# 科学问题：如何同时呈现效应与显著性并突出通路成员或目标基因？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "hyperbolic-volcano-highlights")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复双曲线目标点/背景层级并去掉重复通路图例，六个目标基因标签可读；输入/派生值均为合成示例，不冒充来源结果。
# 顶刊复现 | 双曲线火山图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1020)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/crispr_screen.csv"))
df$product <- abs(df$epsilon)*(-log10(df$pvalue))
## 3. 教程的NTC比例扫描：仅为演示统计量，不宣称校准的FDR ====
cutoffs <- seq(0,max(df$product),by=.1)
scan <- data.frame(threshold=cutoffs,ntc_fraction=NA_real_,hits=NA_integer_)
for (j in seq_along(cutoffs)) {
  hit <- subset(df,product>cutoffs[j]);scan$hits[j] <- nrow(hit)
  scan$ntc_fraction[j] <- sum(hit$NTC)/nrow(hit)
}
thres <- scan$threshold[which(scan$ntc_fraction<.05 & scan$hits>0)[1]]
write.csv(scan,file.path(root,"result/tutorial_ntc_threshold_scan.csv"),row.names=FALSE)
df$base_group <- case_when(df$NTC~"Non-targeting",df$product>thres & df$epsilon>0~"Positive hits",df$product>thres & df$epsilon<0~"Negative hits",TRUE~"Other genes")
df$final <- ifelse(df$Pathway=="",df$base_group,df$Pathway)
cols <- c("Non-targeting"="black","Other genes"="grey65","Positive hits"="thistle","Negative hits"="lightblue","GPI anchor"="#0000cd","Ubiquitin/Proteasome"="#7a15cb","Autophagy"="#228623","Mitochondria"="#e4390d","UFMylation"="#cc7703")
curve <- data.frame(x=c(seq(-1.2,-.015,length.out=400),seq(.015,1.2,length.out=400)),side=rep(c("negative","positive"),each=400));curve$y<-thres/abs(curve$x)
p <- ggplot(df,aes(epsilon,-log10(pvalue))) + geom_point(aes(colour=base_group),size=.85,alpha=.5) +
  geom_point(data=subset(df,Pathway!=""),aes(fill=Pathway),shape=21,colour="black",size=1.8,stroke=.3) +
  geom_line(data=curve,aes(x,y,group=side),inherit.aes=FALSE,linetype=3,linewidth=.35) +
  geom_text_repel(data=subset(df,Label==1),aes(label=gene),size=2.4,seed=60,max.overlaps=Inf,min.segment.length=0) +
  scale_colour_manual(values=cols,breaks=c("Other genes","Non-targeting","Positive hits","Negative hits"),name="Significant") + scale_fill_manual(values=cols,name="Pathway")+guides(colour=guide_legend(order=1),fill=guide_legend(order=2)) + coord_cartesian(xlim=c(-1.2,1.2),ylim=c(0,4.5)) +
  labs(x="Knockdown phenotype (AU)",y="-log10 P",subtitle=sprintf("Tutorial NTC-fraction threshold = %.1f",thres)) + theme(legend.text=element_text(size=6),legend.key.height=unit(3,"mm"),plot.subtitle=element_text(size=7))
save_plot(p,file.path(fig_dir,"20a_hyperbolic_crispr"),160,110)
## 4. 固定乘积阈值不等于RNA-seq的FDR校正 ====
rna <- read.csv(file.path(module_root,"data/rna_volcano.csv"))
rna$group <- case_when(abs(rna$logFC)*(-log10(rna$padj))>6 & rna$logFC>0~"Up",abs(rna$logFC)*(-log10(rna$padj))>6 & rna$logFC<0~"Down",TRUE~"Other")
cc <- data.frame(x=c(seq(-4,-.1,length.out=300),seq(.1,4,length.out=300)),side=rep(1:2,each=300));cc$y<-6/abs(cc$x)
p2 <- ggplot(rna,aes(logFC,-log10(padj))) + geom_point(aes(colour=group),size=.65,alpha=.65) + geom_line(data=cc,aes(x,y,group=side),inherit.aes=FALSE,linetype=3,linewidth=.3) +
  scale_colour_manual(values=c(Up="thistle",Down="lightblue",Other="grey55")) + coord_cartesian(ylim=c(0,15)) + labs(x="log2FC",y="-log10 adjusted P",subtitle="Fixed product cutoff = 6 (illustrative)") + theme(legend.position="none",plot.subtitle=element_text(size=7))
save_plot(p2,file.path(fig_dir,"20b_fixed_hyperbolic_rna"),110,100)
## 5. 常规阈值配合六个目标基因标签 ====
dd <- read.csv(file.path(module_root,"data/conventional.csv"))
dd$group <- case_when(dd$interest~"interest",dd$padj<.05 & abs(dd$logFC)>1~"sig",TRUE~"not sig")
p3 <- ggplot(dd,aes(logFC,-log10(padj))) + geom_point(aes(fill=group,size=group),shape=21,colour="#343434",stroke=.2) +
  scale_fill_manual(values=c(interest="#7a0478",sig="#e2c5d2","not sig"="#d9d9d9")) + scale_size_manual(values=c(interest=2.5,sig=1.5,"not sig"=.7)) +
  geom_vline(xintercept=c(-1,1),linetype=3,linewidth=.3) + geom_hline(yintercept=-log10(.05),linetype=3,linewidth=.3) + geom_label_repel(data=subset(dd,interest),aes(label=gene),size=2.6,fontface="italic",seed=61,max.overlaps=Inf,min.segment.length=0) +
  labs(x="log2FC",y="-log10 adjusted P") + theme(legend.position="none")
save_plot(p3,file.path(fig_dir,"20c_conventional_interest_genes"),110,100)
write.csv(df,file.path(root,"result/crispr_classified.csv"),row.names=FALSE)

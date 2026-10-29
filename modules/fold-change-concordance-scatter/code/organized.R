# Science式FC–FC差异基因散点图 / Science-style fold-change concordance scatter
# 科学问题：两个扰动比较中的表达变化方向是否一致，哪些基因最突出？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "fold-change-concordance-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复中心连续点云、上下调目标符号及标签尺度，保留双FDR编码；输入/派生值均为合成示例，不冒充来源结果。
# 让我看看谁还不会画FC-FC散点图 | Science复现
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1002)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/fold_changes.csv"))
# 仅调整拟态绘图输入：中心黑点云连续，不把DE点集中成两个人工孤岛。
bg<-df$FDR1>.02 & df$FDR2>.02
df$x[bg]<-.65*df$x[bg];df$y[bg]<-.65*df$y[bg]
fdr_cutoff <- .01
top_n <- 100
down <- unique(c(head(df$gene[order(ifelse(df$FDR1<fdr_cutoff,df$x,Inf))],top_n),head(df$gene[order(ifelse(df$FDR2<fdr_cutoff,df$y,Inf))],top_n)))
up <- unique(c(head(df$gene[order(ifelse(df$FDR1<fdr_cutoff,-df$x,Inf))],top_n),head(df$gene[order(ifelse(df$FDR2<fdr_cutoff,-df$y,Inf))],top_n)))
df$regulation <- case_when(df$gene %in% union(down,up) & df$x*df$y<0 ~ "discordant",df$gene %in% down ~ "down",df$gene %in% up ~ "up",TRUE~"n.s.")
up_names<-c("TRIM22","VTCN1","ACE2","LDLRAD1","EDN2","SH2D3C","SLC2A4","ATP6V1B1","CLIC3","EXOC3L","TP53INP1","CHPF","EMP2","BMF")
down_names<-c("KRAS","LSMEM2","HMGA2","ETV1","NT5E","FOSL1","SPRY4","PTX3","TNS4","ETV4","ETV5","DUSP6","EGR1","GPR3")
iu<-head(which(df$regulation=="up")[order(-df$y[df$regulation=="up"])],length(up_names))
id<-head(which(df$regulation=="down")[order(df$y[df$regulation=="down"])],length(down_names))
df$gene[iu]<-up_names[seq_along(iu)];df$gene[id]<-down_names[seq_along(id)]
labels <- df[union(iu,id),]
## 3. FC–FC：点大小和透明度分别编码两个FDR ====
p <- ggplot(df,aes(x,y)) + geom_hline(yintercept=0,colour="grey70",linewidth=.3) + geom_vline(xintercept=0,colour="grey70",linewidth=.3) +
  geom_point(data=subset(df,regulation=="n.s."),size=.48,alpha=.25,shape=21,stroke=.15) +
  geom_point(data=subset(df,regulation!="n.s."),aes(colour=regulation,size=-log10(FDR1),alpha=-log10(FDR2))) +
  scale_colour_manual(values=c(up="#cb181d",down="#2171b5",discordant="#c6a58b")) + scale_size_continuous(range=c(.6,2.2)) + scale_alpha_continuous(range=c(.2,.9)) +
  geom_text_repel(data=labels,aes(label=gene),size=2.65,fontface="italic",max.overlaps=Inf,min.segment.length=0,seed=42,box.padding=.2) +
  coord_cartesian(xlim=c(-3.5,3.5),ylim=c(-6,5)) + labs(x="KRAS vs NS siRNA (log2FC)",y="ERKi 24 h vs Control (log2FC)") +
  theme_bw(base_size=9,base_family="Arial") + theme(panel.grid.major=element_line(colour="grey94",linewidth=.2),panel.grid.minor=element_blank(),legend.position="none")
save_plot(p,file.path(fig_dir,"02_science_fc_fc"),115,110)
write.csv(df,file.path(root,"result/classified_points.csv"),row.names=FALSE)

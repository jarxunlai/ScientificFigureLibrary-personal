# Cell式双比较显著性散点图 / Cell-style dual-comparison significance scatter
# 科学问题：两种比较中共同显著与仅单独显著的基因分布有何不同？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "dual-comparison-significance-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：补回缺失的core-signature完整图，修复连续点云并区分两种比较场景；输入/派生值均为合成示例，不冒充来源结果。
# 让我看看谁还不会画FC-FC散点图 | Cell 复现
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1003)
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
df$regulation <- case_when(df$x>0 & df$y>0 & df$FDR1<fdr_cutoff & df$FDR2<fdr_cutoff~"DE up (both)",df$x<0 & df$y<0 & df$FDR1<fdr_cutoff & df$FDR2<fdr_cutoff~"DE down (both)",df$x>0 & df$y>0 & pmin(df$FDR1,df$FDR2)<fdr_cutoff~"DE up (one)",df$x<0 & df$y<0 & pmin(df$FDR1,df$FDR2)<fdr_cutoff~"DE down (one)",TRUE~"not significant")
labels <- df %>% filter(grepl("both",regulation)) %>% arrange(desc(abs(x)+abs(y))) %>% head(20)
cols <- c("DE up (both)"="#b2182b","DE up (one)"="#efb7bc","DE down (both)"="#21618c","DE down (one)"="#8bbdd5","not significant"="white")
## 3. 区分两个比较中同时/单独显著的基因 ====
p <- ggplot(df,aes(x,y)) + geom_abline(slope=1,intercept=0,linetype=3,linewidth=.3,colour="grey65") +
  geom_point(aes(fill=regulation),shape=21,colour="grey35",stroke=.18,size=.85) +
  geom_hline(yintercept=0,colour="grey50",linewidth=.35) + geom_vline(xintercept=0,colour="grey50",linewidth=.35) +
  geom_text_repel(data=labels,aes(label=gene),size=2.55,max.overlaps=Inf,min.segment.length=0,seed=43,fontface="italic") +
  scale_fill_manual(values=cols) + labs(title="monocytes: severe vs moderate",x="log2FC moderate: Dexa vs ctrl.",y="log2FC severe: Dexa vs ctrl.") +
  coord_cartesian(xlim=c(-3.5,3.5),ylim=c(-5.5,5.5)) + theme_bw(base_size=9,base_family="Arial") + theme(panel.grid=element_blank(),legend.position="inside",legend.position.inside=c(.99,.02),legend.justification=c(1,0),legend.text=element_text(size=6))
# 首个来源图是core-signature视图，不能只输出第二个severe/moderate视图。
core<-df;core$x<-core$x*.65;core$y<-core$y*.42
# 标注按颜色语义选择：来源上调符号只放在红色组，不把IL1R2误标到蓝色下调点。
up_core<-subset(core,regulation=="DE up (both)");down_core<-subset(core,regulation=="DE down (both)")
up_core<-up_core[order(-up_core$y),][1:8,];down_core<-down_core[order(down_core$y),][1:6,]
up_core$gene<-c("IL1R2","TSC22D3","CD163","DDIT4","ZFP36L2","SAP30","JDP2","PER1")
down_core$gene<-c("IL1B","CCL3L3","CCL3","CXCL8","G0S2","CD300E")
core_labels<-rbind(up_core,down_core)
write.csv(core_labels,file.path(base,"result/03_core_visual_labels.csv"),row.names=FALSE)
p_core<-ggplot(core,aes(x,y))+geom_point(aes(fill=regulation),shape=21,colour="grey35",stroke=.18,size=1.05)+
 geom_hline(yintercept=0,linewidth=.35)+geom_vline(xintercept=0,linewidth=.35)+
 geom_text_repel(data=core_labels,aes(label=gene),size=2.7,seed=43,max.overlaps=Inf)+
 geom_abline(slope=1,intercept=0,linetype=2,colour="grey75",linewidth=.35)+
 scale_fill_manual(values=cols,breaks=names(cols))+coord_cartesian(xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))+
 labs(title="monocytes: Dexa vs. ctrl. core signature",x="log2FC moderate: Dexa vs ctrl.",y="log2FC severe: Dexa vs ctrl.")+
 theme_bw(base_size=10)+theme(panel.grid=element_blank(),legend.position="inside",legend.position.inside=c(.99,.02),legend.justification=c(1,0),legend.title=element_blank(),legend.text=element_text(size=5.8))+
 annotate("label",x=.5,y=2.42,label="Dexa core up",fill="#efb7bc",label.size=.2,size=2.4)+
 annotate("label",x=-1.65,y=-2.40,label="Dexa core down",fill="#8bbdd5",label.size=.2,size=2.4)
save_plot(p_core,file.path(fig_dir,"03_core_signature"),115,110)
write.csv(core,file.path(base,"result/03_core_visual_input.csv"),row.names=FALSE)
save_plot(p,file.path(fig_dir,"03_cell_fc_fc"),115,110)
write.csv(df,file.path(root,"result/classified_points.csv"),row.names=FALSE)

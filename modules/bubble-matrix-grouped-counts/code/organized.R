# 突触motif气泡矩阵与分组免疫计数 / Synaptic-motif bubble matrix and grouped immune counts
# 科学问题：类别组成比例和多组细胞数量分别如何清晰对照？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "bubble-matrix-grouped-counts")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复motif气泡面积档位、顶部类别标签，并补回分组计数的均值横线；输入/派生值均为合成示例，不冒充来源结果。
# 为什么你的图总是显得“土气”，而顶刊的图却看起来那么“高级”？
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1014)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/motif_matrix.csv"))
types <- unique(df$M_type)
cols <- setNames(c("#73943c","#5fb946","#75c492","#a8862e","#f3d458","#00679a","#4890cd","#b3c5e4","#5f4a9e","#9f90c4","#07733a","#e71f8e","#8b288e","#c680b6","#7a1131","#d21f27","#f37a37","#8d55a1"),types)
df$M_type <- factor(df$M_type,levels=types);df$Motif <- factor(df$Motif,levels=18:1)
df$size_group <- cut(df$frac_syn,breaks=c(-Inf,.05,.15,.30,Inf),labels=c("<5%","5–15%","15–30%",">=30%"),right=FALSE)
## 3. 同时编码突触比例、细胞形态类别和抑制性motif ====
p1 <- ggplot(subset(df,frac_syn>=.05),aes(M_type,Motif)) + geom_vline(xintercept=c(3.5,5.5,8.5,12.5,14.5,17.5),linewidth=.25,colour="grey55") +
  geom_point(aes(size=size_group,fill=M_type),shape=21,colour="black",stroke=.2) + scale_fill_manual(values=cols,guide="none") +
  scale_size_manual(values=c("<5%"=0,"5–15%"=1.7,"15–30%"=2.7,">=30%"=3.9),name="Frac Syn",drop=FALSE) + scale_x_discrete(position="top") +
  labs(x="Output fraction",y="Inhibitory motif group") + theme(axis.text.x=element_text(angle=90,hjust=0,size=7.8),axis.text.y=element_text(size=7),axis.line=element_blank(),axis.ticks=element_blank())
save_plot(p1,file.path(fig_dir,"14a_motif_bubble_matrix"),150,128)
## 4. 分组细胞计数；空心/实心区别处理而不是重复图例 ====
dd <- read.csv(file.path(module_root,"data/immune_counts.csv"))
dd$Cell <- factor(dd$Cell,levels=unique(dd$Cell));dd$Group <- factor(dd$Group,levels=c("GF mock","SPF mock","GF PCV13","SPF PCV13"))
cols2 <- c("GF mock"="#3853a2","SPF mock"="#31b24c","GF PCV13"="#3853a2","SPF PCV13"="#31b24c")
p2 <- ggplot(dd,aes(Cell,Value,colour=Group,shape=Group,group=Group)) + geom_vline(xintercept=seq(1.5,7.5),linewidth=.3,colour="grey50") +
  stat_summary(fun=mean,geom="crossbar",width=.55,linewidth=.3,position=position_dodge(.7),show.legend=FALSE)+
  geom_point(position=position_jitterdodge(jitter.width=.08,dodge.width=.7,seed=54),size=1.6,stroke=.6) +
  scale_colour_manual(values=cols2) + scale_shape_manual(values=c(1,1,16,16)) + scale_y_continuous(limits=c(0,6000),labels=scales::comma,breaks=seq(0,6000,2000)) +
  labs(x=NULL,y="Cells per sample") + theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="right",legend.text=element_text(size=7))
save_plot(p2,file.path(fig_dir,"14b_grouped_immune_counts"),160,90)

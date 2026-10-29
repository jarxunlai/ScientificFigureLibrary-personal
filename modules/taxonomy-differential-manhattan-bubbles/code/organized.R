# 差异OTU分类曼哈顿气泡图 / Taxonomy-aware differential-OTU Manhattan bubbles
# 科学问题：不同微生物类别在三个比较中的丰度和倍数变化如何分布？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "taxonomy-differential-manhattan-bubbles")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复显著/非显著空实心气泡层级，三比较共享底部分类轴与分类色条；输入/派生值均为合成示例，不冒充来源结果。
# 顶刊复现 | 差异 OTU 曼哈顿图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1017)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/otu.csv"))
orders <- unique(df$Order);families <- unique(df$Family)
cols <- setNames(c("#dd504e","#eb784c","#f3a13d","#f5d46b","#d8f090","#b0d6a7","#69c477","#8bb2c3","#3881ae","#dacce3","#ac87bd","#ebc5de","#ce73b0","#ba4093","#ccdeb6","#a1c273","#658b3c","#f2e7c9","#dec575","#947a23","#87d1d0","#41aba9","#b5c4d7","#879ebe","#9a9192"),orders)
df$Family <- factor(df$Family,levels=families);df$Order <- factor(df$Order,levels=orders)
df$fill <- ifelse(df$significant,cols[df$Order],"white")
panels <- list()
## 3. 三个比较共享分类顺序，点面积映射丰度，实心标记显著 ====
for (contrast in unique(df$Comparison)) {
  dd <- subset(df,Comparison==contrast)
  p <- ggplot(dd,aes(Family,logFC)) + geom_hline(yintercept=0,colour="grey70",linewidth=.25) +
    geom_point(aes(size=Abundance,colour=Order,fill=fill),shape=21,stroke=.3,position=position_jitter(width=.15,height=0,seed=57)) +
    scale_colour_manual(values=cols) + scale_fill_identity() + scale_size_continuous(range=c(.5,4.8),breaks=c(.01,.05,.1)) +
    coord_cartesian(ylim=c(-4,4)) + labs(y="log2FC",x=NULL,title=contrast) + theme_bw(base_size=7) + theme(panel.grid=element_line(colour="grey90",linewidth=.2),axis.text.x=element_text(angle=90,hjust=1,size=5.8),legend.position="none",plot.title=element_text(hjust=1,size=7))
  panels[[contrast]] <- p;save_plot(p,file.path(fig_dir,paste0("17_",gsub("[^A-Za-z0-9]","_",contrast))),180,60)
}
family_keys <- df %>% distinct(Family,Order)
strip <- ggplot(family_keys,aes(Family,1,fill=Order)) + geom_tile(height=.1,width=.9) + scale_fill_manual(values=cols) +
  scale_y_continuous(limits=c(.8,1.2)) + theme_void() + theme(legend.position="none")
ab <- data.frame(x=1:4,y=1,Abundance=c(.001,.01,.05,.1))
abp <- ggplot(ab,aes(x,y,size=Abundance)) + geom_point(shape=21,fill="white",colour="grey30",stroke=.25) +
  geom_text(aes(y=.6,label=Abundance),size=2) + scale_size_continuous(range=c(.5,4.8),limits=c(.001,.1),guide="none") +
  coord_cartesian(ylim=c(.4,1.5)) + labs(title="Relative abundance") + theme_void() + theme(plot.title=element_text(size=6,hjust=.5))
save_plot(abp,file.path(fig_dir,"17_abundance_legend"),65,25)
# 组合图只在最下一个比较显示共享分类轴，避免重复三次标签挤掉主面板。
for(k in 1:2)panels[[k]]<-panels[[k]]+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
combo <- wrap_plots(c(panels,list(strip,abp)),ncol=1,heights=c(1,1,1,.12,.22))
save_plot(combo,file.path(fig_dir,"17_otu_manhattan"),180,165)

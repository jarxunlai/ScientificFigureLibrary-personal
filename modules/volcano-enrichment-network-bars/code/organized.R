# 火山图与富集网络及条形图组合 / Volcano plots with enrichment networks and bars
# 科学问题：基因层面的变化与预计算的通路结果如何联合展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "volcano-enrichment-network-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：修正富集网络颜色强度方向与节点层级，恢复窄火山/宽网络及三图组合比例；输入/派生值均为合成示例，不冒充来源结果。
# 顶刊复现 | 火山图+富集网络图/条形图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1016)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/volcano.csv"));en <- read.csv(file.path(module_root,"data/enrichment.csv"))
df$group <- case_when(df$logFC>.3 & df$pvalue<.05~"Up",df$logFC < -0.3 & df$pvalue<.05~"Down",TRUE~"Not Sig")
nodes <- read.csv(file.path(module_root,"data/network_nodes.csv"));edges <- read.csv(file.path(module_root,"data/network_edges.csv"))
## 3. 案例一：横置火山图 + 预计算富集网络 ====
p1 <- ggplot(df,aes(-log10(pvalue),logFC,colour=group)) + geom_point(size=.45,alpha=.65) + scale_colour_manual(values=c(Up="#c64b49",Down="#2586b8","Not Sig"="grey75")) +
  labs(title="Endothelial cells",x="-log10 FDR",y="log2FC") + theme(legend.position="none")
pn <- ggplot() + ggforce::geom_mark_ellipse(data=nodes,aes(x,y,group=cluster),fill="#f2f2f9",colour="#d5d5e4",alpha=.4,linewidth=.25,expand=unit(1.5,"mm")) +
  geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),colour="#9494d5",alpha=.45,linewidth=.25) +
  geom_point(data=nodes,aes(x,y,size=Size,colour=score),alpha=.9) + scale_colour_gradient(low="#c1c1fa",high="#1717e8",name="-log10 FDR") + scale_size_continuous(range=c(1,3.2),name="Pathway size") +
  geom_label_repel(data=subset(nodes,label==1),aes(x,y,label=Term),size=2.4,linewidth=.2,label.padding=unit(.5,"mm"),seed=56,max.overlaps=Inf) + labs(title="Enrichment network") + theme_void(base_size=7) + theme(legend.position="right",legend.key.height=unit(3,"mm"))
save_plot(p1,file.path(fig_dir,"16a_rotated_volcano"),65,100);save_plot(pn,file.path(fig_dir,"16a_network"),115,100)
save_plot(p1+pn+plot_layout(widths=c(.85,1.85)),file.path(fig_dir,"16a_volcano_network"),210,115)
## 4. 案例二：火山图与上下调富集条形图 ====
pv <- ggplot(df,aes(logFC,-log10(pvalue),colour=group)) + geom_point(size=.55,alpha=.7) + geom_hline(yintercept=-log10(.05),linetype=3,linewidth=.3) + geom_vline(xintercept=c(-.3,.3),linetype=3,linewidth=.3) +
  scale_colour_manual(values=c(Down="#513965",Up="#9e3f58","Not Sig"="grey75")) + labs(title="Plasma proteins expression",x="log2FC (Young vs Aged)",y="-log10 P") + theme(legend.position="none")
bars <- list()
for (g in c("Down","Up")) {
  dd <- subset(en,Direction==g);dd$Term <- factor(dd$Term,levels=rev(dd$Term))
  p <- ggplot(dd,aes(score,Term)) + geom_col(fill=c(Down="#665077",Up="#b87a8b")[g],width=.58) + geom_text(aes(x=.1,label=Term),hjust=0,size=2.3) +
    labs(x="-log10 P",y=NULL,title=paste(g,"proteins")) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=8))
  bars[[g]] <- p;save_plot(p,file.path(fig_dir,paste0("16b_",g,"_enrichment")),70,100)
}
save_plot(pv,file.path(fig_dir,"16b_volcano"),85,100)
save_plot(pv+bars$Down+bars$Up+plot_layout(widths=c(1.15,1,1)),file.path(fig_dir,"16b_volcano_bars"),215,100)
## 5. 案例三：带标签的火山图 + 正负方向富集条形图 ====
pv3 <- pv + scale_colour_manual(values=c(Down="#66cc00",Up="#ffcc00","Not Sig"="grey75")) + labs(title="Differential plasma proteins") +
  geom_text_repel(data=df %>% arrange(pvalue) %>% head(16),aes(label=gene),size=2.2,seed=57,max.overlaps=Inf)
en$signed <- ifelse(en$Direction=="Up",en$score,-en$score);en$Term <- factor(en$Term,levels=rev(en$Term))
pb3 <- ggplot(en,aes(signed,Term,fill=Direction)) + geom_col(width=.7,colour="grey35",linewidth=.15) + geom_vline(xintercept=0,linewidth=.3) + scale_fill_manual(values=c(Down="#66cc00",Up="#ffcc00")) +
  labs(x="signed -log10 FDR",y=NULL) + theme(axis.text.y=element_text(size=5),legend.position="bottom")
save_plot(pv3,file.path(fig_dir,"16c_labeled_volcano"),90,110);save_plot(pb3,file.path(fig_dir,"16c_divergent_bars"),115,110)
save_plot(pv3+pb3+plot_layout(widths=c(1,1.25)),file.path(fig_dir,"16c_volcano_divergent"),180,110)

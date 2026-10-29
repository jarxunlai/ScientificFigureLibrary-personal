# 三类多分组火山图 / Three multi-group volcano layouts
# 科学问题：多个细胞群或比较的差异表达模式和极值基因如何同图比较？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "multigroup-volcano-layouts")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 三类多分组火山图：按来源保留整体+九宫格、横向九组、中心色条三幅独立图。
# 原始 CSV 仍为合成示例；拟态效应量和P值另存 result/，不冒充论文结果。
# 使用 source/author-code-extracted.txt 核对布局/颜色，未执行安装或上游分析代码。
Sys.setlocale("LC_CTYPE","English_United States.utf8")



source(file.path(module_root,"code/export.R"))
options(warn=1)
options(device=function(...) grDevices::cairo_pdf(file=tempfile(fileext=".pdf"),family="Arial",...))
library(ggrepel)
fig_dir <- file.path(base,"figure"); res <- file.path(base,"result")
set.seed(260115)
clusters <- c("B","CD14+ Mono","CD8 T","DC","FCGR3A+ Mono","Memory CD4 T","Naive CD4 T","NK","Platelet")
cols <- setNames(c("#76a7c2","#80b97d","#795fa3","#edac77","#5f91ac","#bc8763","#a1a452","#d9c27e","#d14a3e"),clusters)
cols2 <- setNames(c("#c0322f","#26aee7","#96ca39","#51aa44","#f7bcd0","#f1bd58","#68c2b2","#f74141","#b783b6"),clusters)
# 文字来自来源图中的基因名；下面所有数值仍是合成，名称不是验证结果。
label_up <- list(c("CD79B","TCL1A","HLA-DQB1","CD74","HLA-DQA1"),
  c("LYZ","TYROBP","CST3","MS4A6A","CD14"),c("CCL5","B2M","NKG7","GZMH","CST7"),
  c("FCER1A","HLA-DPB1","HLA-DQB1","HLA-DQA1","HLA-DPA1"),
  c("LST1","FCGR3A","FCER1G","AIF1","MS4A7"),c("CD2","AQP3","LTB","IL7R","IL32"),
  c("CCR7","LEF1","LDHB","NOSIP","PRKCQ-AS1"),c("GNLY","GZMB","FGFBP2","NKG7","PRF1"),
  c("PPBP","PF4","GNG11","CLU","SDPR"))
label_down <- list(c("NKG7","TYROBP","S100A9","S100A8","LYZ"),c("CCL5","PTPRCAP","IL32","NKG7","LTB"),
  c("HLA-DRA","CST3","LYZ","S100A8","S100A9"),c("CD3E","CD7","LTB","IL32","CD3D"),
  c("LTB","CD3D","CCL5","IL32","PTPRCAP"),c("LYZ","TYROBP","S100A8","S100A9","HLA-DRA"),
  c("LYZ","TYROBP","CST3","S100A8","HLA-DRA"),c("LYZ","CST3","HLA-DRA","S100A8","S100A9"),
  c("RPL14","MALAT1","RPL9","RPS24","RPS4Y1"))
raw <- read.csv(file.path(module_root,"data/group_volcano.csv"))
raw$cluster <- factor(raw$cluster,levels=clusters)
theme_source <- theme_classic(base_size=10,base_family="sans")+
  theme(axis.text=element_text(colour="black"),axis.line=element_line(linewidth=.35),
        axis.ticks=element_line(linewidth=.35),legend.title=element_blank(),plot.margin=margin(3,3,3,3))

## 1. 整体+九宫格：恢复组别色条、紧邻的固定坐标、中央白带 ====
a <- list()
for(k in seq_along(clusters)) {
  dd <- raw[raw$cluster==clusters[k],]
  take <- c(45,27,30,52,22,29,26,36,86)[k]
  dd <- dd[seq_len(take),]
  u <- (rank(dd$logFC)-.5)/nrow(dd)
  direction <- ifelse(u>.5,1,-1)
  center <- c(.65,.6,.7,.95,.65,.7,.6,1.2,2.6)[k]
  dd$logFC <- direction*(.51+center*abs(qnorm(pmin(.98,pmax(.02,u))))+.1*runif(nrow(dd)))
  dd$padj <- 10^(-runif(nrow(dd),1.3,4.8))
  a[[k]] <- dd
}
first <- bind_rows(a); first$cluster <- factor(first$cluster,levels=clusters)
write.csv(first,file.path(res,"15a_visual_example.csv"),row.names=FALSE)
background <- list(annotate("rect",xmin=-4,xmax=-.5,ymin=-Inf,ymax=Inf,fill="#C7DFE8"),
                   annotate("rect",xmin=.5,xmax=4,ymin=-Inf,ymax=Inf,fill="#FFD8DF"))
p <- ggplot(first,aes(logFC,-log10(padj)))+background+
  geom_point(aes(colour=cluster),size=1.25,alpha=.8)+scale_colour_manual(values=cols)+
  scale_x_continuous(breaks=c(-4,-2,-.5,.5,2,4),expand=c(0,0))+
  scale_y_continuous(breaks=c(1.3,2.3,3.3,4.3),expand=c(0,0))+
  coord_cartesian(xlim=c(-4.4,4.4),ylim=c(1.15,5.05))+
  labs(x=expression(log[2](FC)),y=expression(-log[10](adj.p.value)))+
  theme_source+theme(legend.position="none")
# 原生 ggplot 分面 strip 逐项着色，不需要 scRNAtoolVis/ggh4x。
pf <- p+facet_wrap(~cluster,ncol=3)+
  scale_x_continuous(breaks=c(-4,4),expand=c(0,0))+
  scale_y_continuous(breaks=c(1.3,2.8,4.3),expand=c(0,0))+
  labs(y=NULL)+theme(panel.spacing=unit(.75,"mm"),strip.text=element_text(colour="white",face="bold",size=8),
    strip.background=element_rect(fill="#A8B1B4",colour=NA),axis.text=element_text(size=7))
g <- ggplotGrob(pf)
strip_ids <- grep("strip-t",g$layout$name)
strip_ids <- strip_ids[order(g$layout$t[strip_ids],g$layout$l[strip_ids])]
for(k in seq_along(strip_ids)) {
  strip <- g$grobs[[strip_ids[k]]]$grobs[[1]]
  bg <- grep("strip.background",names(strip$children))
  strip$children[[bg]]$gp$fill <- scales::alpha(cols[k],.72)
  g$grobs[[strip_ids[k]]]$grobs[[1]] <- strip
}
save_plot(p,file.path(fig_dir,"15a_overlaid_volcano"),112,80)
save_plot(g,file.path(fig_dir,"15a_nine_facets"),118,80)
combined_a <- cowplot::plot_grid(p,g,nrow=1,rel_widths=c(1.12,1),align="hv",axis="tb")
save_plot(combined_a,file.path(fig_dir,"15a_overview_and_facets"),240,86)

## 2. 横向九组：来源中横轴到300、纵轴效应到9；不再只画灰色散点云 ====
b <- list(); tops <- list()
for(k in seq_along(clusters)) {
  dd <- raw[raw$cluster==clusters[k],]
  u <- (rank(dd$logFC)-.5)/nrow(dd)
  z <- qnorm(u)
  up_limit <- c(4,4.7,3.0,4.2,3.7,2.0,2.1,5.3,8.6)[k]
  down_limit <- c(4.1,3.7,4.0,3.5,3.6,4.0,4.2,4.1,4.2)[k]
  dd$logFC <- ifelse(z>0,up_limit*(z/max(z))^.75,-down_limit*(abs(z)/abs(min(z)))^.75)
  logp <- (abs(dd$logFC)+.06)^1.65*runif(nrow(dd),.68,1.22)
  logp <- logp/max(logp)*c(275,280,240,265,245,240,230,280,300)[k]
  dd$padj <- 10^(-logp)
  dd$gene <- paste0("Example_",k,"_",seq_len(nrow(dd)))
  hi <- order(dd$logFC,decreasing=TRUE)[1:5]; lo <- order(dd$logFC)[1:5]
  dd$gene[hi] <- label_up[[k]]; dd$gene[lo] <- label_down[[k]]
  tops[[k]] <- dd[c(hi,lo),]; b[[k]] <- dd
}
second <- bind_rows(b); top <- bind_rows(tops)
write.csv(second,file.path(res,"15b_visual_example.csv"),row.names=FALSE)
p2 <- ggplot(second,aes(-log10(padj),logFC))+
  geom_point(colour="#BDBDBD",size=.47,alpha=.9)+
  geom_hline(yintercept=c(-.05,.05),linetype="dashed",linewidth=.25,colour="#888888")+
  geom_point(data=top,aes(colour=cluster),size=1.1)+
  geom_text_repel(data=top,aes(label=gene,colour=cluster),size=2.05,fontface="italic",
    seed=233,max.overlaps=Inf,min.segment.length=0,box.padding=.15,
    force=4,force_pull=1,segment.alpha=.6,segment.size=.2,max.time=2)+
  facet_grid(.~cluster)+scale_colour_manual(values=cols2)+
  scale_x_continuous(breaks=c(0,100,200,300),limits=c(0,330),expand=c(0,0))+
  scale_y_continuous(breaks=seq(-4,8,2),limits=c(-4.7,9.3),expand=c(0,0))+
  labs(x=expression(-log[10](p_val_adj)),y="avg_log2FC")+
  theme_source+theme(legend.position="none",panel.border=element_rect(fill=NA,linewidth=.3),
    strip.background=element_rect(fill="#D8D8D8",colour="#666666",linewidth=.25),
    strip.text=element_text(size=8,face="bold"),panel.spacing=unit(1.25,"mm"),
    axis.line=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=7))
save_plot(p2,file.path(fig_dir,"15b_horizontal_facets"),270,104)

## 3. 临床分组：恢复组名在色条中、每组灰背景、上下调两侧和10个标签 ====
raw3 <- read.csv(file.path(module_root,"data/group_volcano3.csv"))
groups <- unique(raw3$cluster)
col3 <- setNames(c("#E89676","#EFA888","#F7B99C","#FCD1B9","#FBE6DA","#F7F7F7","#D8E9F1","#85BDD4"),groups)
cases <- list(); top3 <- list()
for(k in seq_along(groups)) {
  dd <- raw3[raw3$cluster==groups[k],]
  u <- (rank(dd$logFC)-.5)/nrow(dd)
  dd$logFC <- ifelse(u<.5,-1.05-abs(qnorm(u))*c(1.0,1.1,1.2,1.0,1.0,1.1,1.0,1.1)[k],
                    1.05+abs(qnorm(u))*c(2.7,1.3,1.3,1.3,1.1,1.0,1.0,1.0)[k])
  dd$y <- k+runif(nrow(dd),-.34,.34)
  dd$gene <- paste0("Example_",seq_len(nrow(dd)))
  hi <- order(dd$logFC,decreasing=TRUE)[1:5];lo <- order(dd$logFC)[1:5]
  dd$gene[hi] <- label_up[[c(9,8,5,4,1,7,6,3)[k]]]
  dd$gene[lo] <- label_down[[c(9,8,5,4,1,7,6,3)[k]]]
  cases[[k]] <- dd;top3[[k]] <- dd[c(hi,lo),]
}
third <- bind_rows(cases);top3 <- bind_rows(top3)
write.csv(third,file.path(res,"15c_visual_example.csv"),row.names=FALSE)
bands <- third %>% group_by(cluster) %>% summarise(y=first(y)%/%1+.5,xmin=min(logFC)-.25,xmax=max(logFC)+.25,.groups="drop")
bands$y <- match(bands$cluster,groups)
center <- data.frame(cluster=groups,y=seq_along(groups),x=0)
p3 <- ggplot(third,aes(logFC,y))+
  geom_rect(data=bands,aes(xmin=xmin,xmax=xmax,ymin=y-.45,ymax=y+.45),inherit.aes=FALSE,fill="#EBEBEB")+
  geom_tile(data=center,aes(x=x,y=y,fill=cluster),inherit.aes=FALSE,width=.5,height=1,colour="#555555",linewidth=.25)+
  geom_point(aes(colour=ifelse(logFC>0,"Significant Up","Significant Down")),size=.72)+
  geom_text_repel(data=top3,aes(label=gene),size=2.7,fontface="italic",seed=233,
    min.segment.length=0,box.padding=.1,force=2,max.overlaps=Inf,segment.size=.3,max.time=3)+
  geom_label(data=center,aes(x=x,y=y,label=cluster),inherit.aes=FALSE,size=2.75,
    label.padding=unit(.4,"mm"),linewidth=.25,fill="white")+
  scale_fill_manual(values=col3,guide="none")+
  scale_colour_manual(values=c("Significant Down"="#0099CC","Significant Up"="#CC3333"))+
  scale_y_continuous(breaks=NULL,expand=c(0,0),limits=c(.45,8.55))+
  scale_x_continuous(breaks=seq(-5,7.5,2.5),expand=expansion(mult=c(.03,.08)))+
  labs(x="Log2 Fold Change",y="Group")+theme_source+
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position="inside",
    legend.position.inside=c(.84,.85),legend.background=element_rect(fill="white",colour=NA),
    legend.text=element_text(size=9),legend.key.size=unit(4,"mm"))+
  guides(colour=guide_legend(override.aes=list(size=2.7)))
save_plot(p3,file.path(fig_dir,"15c_multigroup_jitter_volcano"),210,153)
cat("VISUAL_REPAIR: 3 source layouts; 5 PNG/PDF pairs; synthetic tables saved in result/\n")

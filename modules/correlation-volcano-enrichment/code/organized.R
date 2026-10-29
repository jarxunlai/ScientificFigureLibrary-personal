# 相关性火山图与运行富集曲线 / Correlation volcanoes and enrichment
# 科学问题：哪些基因或蛋白与连续表型相关，一组通路成员是否集中在相关排序的一端？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "correlation-volcano-enrichment")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复纵向火山-运行富集-rug三层比例，增大目标点与标签但保留实际富集计算；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. Pearson相关火山图：P值来自r与样本量的双尾t检验 ====
d <- read.csv(file.path(dat,"wnt_correlations.csv"))
d <- d[order(d$PearsonR),];d$logP <- -log10(d$PValue)
q_cut <- .05; r_cut <- -.3
p_cut <- max(d$PValue[d$FDR<q_cut]); y_cut <- -log10(p_cut)
lab <- subset(d,Gene %in% c("MESD","WLS","TCF7L2","FZD5","LGR4","CTNNB1","BCL9L","CSNK2B","WNT10B"))
lab$segment <- ifelse(lab$FDR<q_cut,"#d66a6a","grey55")
p1 <- ggplot(d,aes(PearsonR,logP))+geom_hline(yintercept=y_cut,linetype=2,colour="#e58a8a",linewidth=.4)+
  geom_point(data=subset(d,!Wnt_Signaling),shape=21,size=1.0,stroke=.15,fill="grey75",colour="grey35",alpha=.7)+
  geom_point(data=subset(d,Wnt_Signaling),aes(fill="WNT signalling"),shape=21,size=1.2,stroke=.2,colour="grey35",alpha=.85)+
  scale_fill_manual(values=c(`WNT signalling`="#07e0e7"),name=NULL)+
  ggrepel::geom_text_repel(data=lab,aes(label=Gene,segment.colour=segment),fontface="italic",size=3.0,max.overlaps=Inf,min.segment.length=0,seed=123,box.padding=.45)+
  scale_x_continuous(limits=c(-.6,.55),breaks=seq(-.4,.4,.2))+labs(x=NULL,y="-log10(P)")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="inside",legend.position.inside=c(.73,.91),plot.margin=margin(3,6,0,5))
## 2. 从全排序列表计算加权运行富集；不把插值曲线冒充GSEA检验 ====
hit <- d$Wnt_Signaling; weight <- abs(d$PearsonR)
step <- ifelse(hit,weight/sum(weight[hit]),-1/sum(!hit)); d$RunningES <- cumsum(step)
write.csv(d,file.path(root,"result","wnt_with_running_es.csv"),row.names=FALSE)
p2 <- ggplot(d,aes(PearsonR,RunningES))+geom_hline(yintercept=0,colour="grey60",linewidth=.3)+geom_line(colour="#00cbd4",linewidth=.45)+
  scale_x_continuous(limits=c(-.6,.55),breaks=seq(-.4,.4,.2))+labs(x=NULL,y="Enrichment\nscore")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=6),axis.title.y=element_text(size=7),plot.margin=margin(0,6,0,5))
p3 <- ggplot(subset(d,hit),aes(PearsonR))+geom_segment(aes(xend=PearsonR,y=0,yend=1),colour="#00cbd4",linewidth=.3)+
  scale_x_continuous(limits=c(-.6,.55),breaks=seq(-.4,.4,.2))+scale_y_continuous(expand=c(0,0))+
  labs(x="Pearson's correlation",y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(0,6,4,5))
out <- p1/p2/p3+plot_layout(heights=c(6.5,1,.45))
save_plot(out,file.path(fig,"14a_pearson_wnt_enrichment"),140,195)
## 3. Kendall tau来自模拟配对蛋白丰度的真正秩相关检验 ====
k <- read.csv(file.path(dat,"kendall_correlations.csv"));k$logP <- -log10(k$PValue)
k$Significance <- ifelse(k$FDR<.05,"Significant","NS")
labs_k <- subset(k,Gene %in% c("HSPA2","LGALS1","MYBPC2","PDLIM3","ARHGAP10","LDHA","NDUFB2","LDHB","BDH1","STEAP4"))
pk <- ggplot(k,aes(tau,logP))+geom_point(data=subset(k,Significance=="NS"),colour="grey75",size=.3,alpha=.25)+
  geom_point(data=subset(k,Significance=="Significant"),colour="#991f1f",size=1,alpha=.4)+
  ggrepel::geom_text_repel(data=labs_k,aes(label=Gene),colour="#991f1f",size=2.8,seed=123,max.overlaps=Inf,box.padding=.5,min.segment.length=0,segment.colour="#991f1f")+
  scale_x_continuous(limits=c(-.55,.55),breaks=c(-.5,-.25,0,.25,.5))+labs(x=expression(tau),y="-log10(P)")
save_plot(pk,file.path(fig,"14b_kendall_tau_volcano"),170,155)

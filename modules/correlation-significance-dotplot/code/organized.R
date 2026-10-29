# 显著性与丰度双编码相关点阵 / Correlation dotplots: P value and abundance
# 科学问题：多组特征间的关联方向、统计显著性和丰度如何分别表达？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "correlation-significance-dotplot")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复密集圆点矩阵比例和黑色描边，区分P值大小与丰度大小图例；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. WGD分组：点大小为-log10P，填色为Spearman rho ====
d <- read.csv(file.path(dat,"wgd_correlations.csv"))
features <- c("MHC-II","MHC-I","Interferon MHC II (2)","Interferon MHC II (1)","Complement","Apoptosis","TNFa via NFKB","Inflam. Resp.","Protein Secretion","IFNa","IFNg (H)","STING1","CIN","Inflam.","IFNg (B)","ISG")
d$Feature <- factor(d$Feature,levels=features)
d$CNV <- factor(d$CNV,levels=c("Chrom. Gain","Chrom. Loss","Arm Gain","Arm Loss","Segment Gain","Segment Loss"))
d$Type <- factor(d$Type,levels=c("Curated","Gene","Hallmark pathways","Malignant metaprograms","MHC"))
d$Group <- factor(d$Group,levels=c("WGD-low","WGD-high"))
d$star <- ifelse(d$pval<.001,"***",ifelse(d$pval<.01,"**",ifelse(d$pval<.05,"*","")))
p <- ggplot(d,aes(CNV,Feature))+
  geom_point(aes(size=-log10(pval),fill=rho),shape=21,colour="black",stroke=.25)+geom_text(aes(label=star),size=2.1,vjust=.70)+
  facet_grid(Type~Group,scales="free_y",space="free_y",switch="y")+
  scale_fill_gradient2(low="#8b81bd",mid="white",high="#822425",limits=c(-.6,.6),name="Spearman rho")+
  scale_size_continuous(range=c(1.2,4.8),name="-log10(P)",breaks=c(1,2,3))+
  labs(x=NULL,y=NULL)+theme_bw(base_size=9,base_family="Arial")+
  theme(panel.grid=element_blank(),strip.background=element_blank(),strip.placement="outside",strip.text.y.left=element_text(angle=0,size=7),axis.text.y=element_text(size=7),legend.text=element_text(size=7),legend.title=element_text(size=8),axis.text.x=element_text(angle=52,hjust=1),panel.spacing.y=unit(1.2,"mm"))
save_plot(p,file.path(fig,"09a_wgd_significance_dotplot"),195,158)
## 2. ARG点阵：紧凑矩阵、三档丰度圆点与中心显著性星号 ====
# 2026-10-27 单图精修：不改变现有合成CSV的分类、rho、P或丰度。
# 来源参考：source/reference-06.png；原图15行×10列，矩阵近方形，圆点应占据足够行高。
# 当前合成分类保留本表名称，不把Proteobacteria冒充为Beta/Gamma等亚类。
a <- read.csv(file.path(dat,"arg_correlations.csv"))
a$dataset <- factor(a$dataset,levels=c("16S-ARG","Shotgun-ARG","16S-Beta-lactam","Shotgun-Beta-lactam","16S-Glycopeptide","Shotgun-Glycopeptide","16S-Tetracycline","Shotgun-Tetracycline","16S-Aminoglycoside","Shotgun-Aminoglycoside"))
ph <- a %>% group_by(Phylum) %>% summarise(ab=mean(abundance),.groups="drop") %>% arrange(ab)
a$Phylum <- factor(a$Phylum,levels=ph$Phylum)
abundance_breaks <- c(.001,.01,.1,1)
abundance_labels <- c("(0.001,0.01]","(0.01,0.1]","(0.1,1]")
abundance_sizes <- setNames(c(3,4,5),abundance_labels)
a$abundance_bin <- cut(a$abundance,breaks=abundance_breaks,labels=abundance_labels,right=TRUE)
a$star <- ifelse(a$pval<.05,"*","")
write.csv(a,file.path(entry,"result/09b_arg_plot_input.csv"),row.names=FALSE)
write.csv(data.frame(Phylum=rev(ph$Phylum),mean_abundance=rev(ph$ab),row_from_top=seq_len(nrow(ph))),
          file.path(entry,"result/09b_arg_row_order.csv"),row.names=FALSE)

p2 <- ggplot(a,aes(dataset,Phylum))+
  geom_vline(xintercept=c(2.5,4.5,6.5,8.5),linewidth=.55,colour="black")+
  geom_point(aes(size=abundance_bin,fill=rho),shape=21,colour="black",stroke=.65)+
  geom_text(aes(label=star),size=3.05,vjust=.78,fontface="bold",show.legend=FALSE)+
  scale_fill_distiller(palette="RdBu",direction=1,limits=c(-.8,.8),
                       breaks=seq(-.8,.8,.4),name="Cor. with ARG")+
  scale_size_manual(values=abundance_sizes,drop=FALSE,name="Relative abundance")+
  scale_x_discrete(expand=expansion(add=.55))+
  scale_y_discrete(expand=expansion(add=.55))+
  coord_fixed(ratio=.65,clip="off")+
  labs(x=NULL,y=NULL)+theme_classic(base_size=10,base_family="Arial")+
  theme(axis.line=element_blank(),panel.border=element_rect(fill=NA,colour="black",linewidth=.65),
        axis.ticks=element_line(colour="grey25",linewidth=.45),axis.ticks.length=unit(1.2,"mm"),
        axis.text=element_text(colour="grey25"),axis.text.y=element_text(size=9.4),
        axis.text.x=element_text(angle=33,hjust=1,vjust=1,size=9.1,margin=margin(t=2)),
        legend.position="right",legend.box="vertical",legend.box.just="left",
        legend.title=element_text(size=10),legend.text=element_text(size=9),
        legend.key=element_rect(fill="white",colour=NA),legend.key.height=unit(6,"mm"),
        legend.key.width=unit(6,"mm"),legend.spacing.y=unit(6,"mm"),
        legend.box.spacing=unit(5,"mm"),plot.margin=margin(4,5,4,3,"mm"))+
  guides(fill=guide_colourbar(order=1,title.position="top",barheight=unit(31,"mm"),
                             barwidth=unit(4,"mm"),ticks.colour="grey60"),
         size=guide_legend(order=2,title.position="top",override.aes=list(fill="white",colour="black",stroke=.65)))
save_plot(p2,file.path(fig,"09b_arg_abundance_dotplot"),178,118)

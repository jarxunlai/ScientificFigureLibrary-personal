# 火山图与顶端内嵌富集条形图 / Volcano plots with top-inset enrichment bars
# 科学问题：上下调基因的变化与其功能条目如何在同一画面关联？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "volcano-inset-enrichment-bars")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：纠正富集条形为火山坐标框内的真正inset，恢复标签和两种阈值版本的共同布局；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 火山图 + 内嵌富集条形图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1024)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/volcano.csv"));en<-read.csv(file.path(module_root,"data/enrichment.csv"))
fc_cutoff <- .58;p_cutoff <- .05;product_cutoff <- 6
bars <- list()
## 3. 预计算富集条形图；不运行GO/KEGG富集分析 ====
for (g in c("Down","Up")) {
  dd<-subset(en,Direction==g);dd$Term<-factor(dd$Term,levels=rev(dd$Term))
  bars[[g]]<-ggplot(dd,aes(score,Term)) + geom_col(fill=c(Down="#e7b2bd",Up="#b9c9d2")[g],width=.55) +
    geom_text(aes(x=.08,label=Term),size=2.2,hjust=0) + geom_text(aes(x=.08,label=Hits),size=1.8,hjust=0,vjust=2.2) +
    labs(x="-log10 P",y=NULL) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=5),axis.title.x=element_text(size=6),plot.margin=margin(1,3,1,3))
  save_plot(bars[[g]],file.path(fig_dir,paste0("24_",g,"_enrichment")),80,60)
}
# 第二个来源是另一组KEGG通路，不可复用第一张GO条形标签与颜色。
kegg_terms<-list(Down=c("GPI-anchor biosynthesis","Basal transcription factors","Terpenoid backbone biosynthesis","Autophagy","Ribosome biogenesis in eukaryotes","Peroxisome","Ribosome","Pantothenate and CoA biosynthesis","RNA transport","Protein export"),
 Up=c("Oxidative phosphorylation","Ubiquitin mediated proteolysis","mRNA surveillance pathway","Autophagy","Spliceosome","mTOR signaling pathway","Porphyrin and chlorophyll metabolism","Other types of O-glycan biosynthesis","Biosynthesis of unsaturated fatty acids","Alzheimer disease"))
kegg_bars<-list();kegg_input<-list()
for(g in names(kegg_terms)){
 seed_scores<-subset(en,Direction==g)$score
 dd<-data.frame(Term=kegg_terms[[g]],Direction=g,score=approx(seq_along(seed_scores),seed_scores,xout=seq(1,length(seed_scores),length.out=10))$y)
 dd$Term<-factor(dd$Term,levels=rev(dd$Term));kegg_input[[g]]<-dd
 kegg_bars[[g]]<-ggplot(dd,aes(score,Term))+geom_col(fill=c(Down="lightblue",Up="thistle")[g],width=.55)+
  geom_text(aes(x=.06,label=Term),hjust=0,size=2.1)+labs(x="-log10 P",y=NULL)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=5),axis.title.x=element_text(size=6),plot.margin=margin(1,3,1,3))
 save_plot(kegg_bars[[g]],file.path(fig_dir,paste0("24_",g,"_kegg_enrichment")),80,65)
}
write.csv(bind_rows(kegg_input),file.path(base,"result/24_kegg_visual_synthetic_input.csv"),row.names=FALSE)
## 4. 传统阈值与双曲线阈值：两种完整版本 ====
for (mode in c("linear","hyperbolic")) {
  if(mode=="linear") {
    df$group<-case_when(abs(df$logFC)>fc_cutoff & df$pvalue<p_cutoff~"Top sides",df$pvalue<p_cutoff~"Top middle",abs(df$logFC)>fc_cutoff~"Bottom sides",TRUE~"Bottom middle")
    cols<-c("Top sides"="#ce0606","Top middle"="#6ca5cc","Bottom sides"="#626262","Bottom middle"="#b1b1b1")
  } else {
    df$group<-case_when(abs(df$logFC)*(-log10(df$pvalue))>product_cutoff & df$logFC>0~"Up",abs(df$logFC)*(-log10(df$pvalue))>product_cutoff & df$logFC<0~"Down",TRUE~"Other")
    cols<-c(Up="thistle",Down="lightblue",Other="#626262")
  }
  p<-ggplot(df,aes(logFC,-log10(pvalue)))
  if(mode=="linear") p<-p + annotate("rect",xmin=-4,xmax=-fc_cutoff,ymin=-log10(p_cutoff),ymax=7,fill="#faf0f2") + annotate("rect",xmin=fc_cutoff,xmax=4,ymin=-log10(p_cutoff),ymax=7,fill="#dde0e4")
  if(mode=="hyperbolic") {
    fill_df<-data.frame(x=seq(product_cutoff/7,4,length.out=300));fill_df$ymin<-product_cutoff/fill_df$x
    p<-p + geom_ribbon(data=fill_df,aes(x=x,ymin=ymin,ymax=7),inherit.aes=FALSE,fill="#f2e5f0") + geom_ribbon(data=fill_df,aes(x=-x,ymin=ymin,ymax=7),inherit.aes=FALSE,fill="#e3eef4")
  }
  p<-p +
    geom_point(aes(colour=group),size=.65,alpha=.7) + scale_colour_manual(values=cols) +
    geom_text_repel(data=df %>% arrange(pvalue) %>% head(16),aes(label=gene),size=2.1,seed=64,max.overlaps=Inf,min.segment.length=0) +
    coord_cartesian(xlim=c(-4,4),ylim=c(0,7)) + labs(x="log2FC",y="-log10 P") + theme(legend.position="none")
  if(mode=="linear") p<-p + geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),linetype=3,linewidth=.3) + geom_hline(yintercept=-log10(p_cutoff),linetype=3,linewidth=.3)
  if(mode=="hyperbolic") {
    cc<-data.frame(x=c(seq(-4,-.1,length.out=300),seq(.1,4,length.out=300)),side=rep(1:2,each=300));cc$y<-product_cutoff/abs(cc$x)
    p<-p+geom_line(data=cc,aes(x,y,group=side),inherit.aes=FALSE,linetype=3,linewidth=.3)
  }
  save_plot(p,file.path(fig_dir,paste0("24_",mode,"_volcano")),130,110)
  # 两个富集条形真正内嵌于同一火山坐标框的上角，不另拼一行冒充inset。
  current_bars<-if(mode=="linear")bars else kegg_bars
  p_full<-p+scale_y_continuous(breaks=c(0,2,4,6))+coord_cartesian(xlim=c(-4,4),ylim=c(0,10),expand=FALSE)
  final<-p_full+inset_element(current_bars$Down,left=.01,bottom=.70,right=.49,top=.99,align_to="panel")+
    inset_element(current_bars$Up,left=.51,bottom=.70,right=.99,top=.99,align_to="panel")
  save_plot(final,file.path(fig_dir,paste0("24_",mode,"_volcano_inset_bars")),145,155)
}

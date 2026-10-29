# 增强子激活抑制双热图 / Enhancer activation and repression dual heatmaps
# 科学问题：转录因子及其组合在不同细胞状态下偏向激活、抑制还是双重作用？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "enhancer-activation-repression-heatmap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复来源细胞状态从上到下顺序和三行共享列对齐，激活/抑制仍各用正确色标；输入/派生值均为合成示例，不冒充来源结果。
# 增强子激活/抑制：两个完整变体 + 三个独立面板。作者代码保留且未执行。
# 输入为模拟序列计数，二项检验仅作可替换示例；显著性以 BH padj<0.05 标记。
library(ggplot2)
library(grid)
library(ggnewscale)


source(file.path(module_root,"code/export.R"))

## 1. 细胞状态变体：数据、顺序、计数检验与摘要 ====
d<-read.csv(file.path(module_root,"data/factor_state_counts.csv"))
d$TFcombo<-factor(d$TFcombo,levels=unique(d$TFcombo))
d$CellState<-factor(d$CellState,levels=c("MegEry","Basophil","Eosinophil","Monocyte","Monocyte P.","Neutrophil","Immature"))
d$Faceting<-factor(d$Faceting,levels=c("Set 1 (single f.)","Set 1 (pairs)","Set 2 (single factors)"))
background_active<-.08; background_repressed<-.03
d$P_active<-p.adjust(pbinom(d$Active-1,d$N,background_active,lower.tail=FALSE),"BH")
d$P_repressed<-p.adjust(pbinom(d$Repressed-1,d$N,background_repressed,lower.tail=FALSE),"BH")
d$PercentActive<-100*d$Active/d$N;d$PercentRepressed<-100*d$Repressed/d$N
summary<-data.frame()
for(tf in levels(d$TFcombo)) {
  block<-d[d$TFcombo==tf,]
  a<-any(block$P_active<.05); r<-any(block$P_repressed<.05)
  class<-if(a&r) "Dual" else if(a) "Activator" else if(r) "Repressor" else "Other"
  summary<-rbind(summary,data.frame(TFcombo=tf,Faceting=as.character(block$Faceting[1]),Class=class))
}
summary$TFcombo<-factor(summary$TFcombo,levels=levels(d$TFcombo))
summary$Faceting<-factor(summary$Faceting,levels=levels(d$Faceting))
summary$Class<-factor(summary$Class,levels=c("Activator","Dual","Repressor","Other"))
write.csv(d,file.path(draft,"result/08a_factor_state_statistics.csv"),row.names=FALSE)
write.csv(summary,file.path(draft,"result/08a_factor_classes.csv"),row.names=FALSE)

## 2. 上下热图与底部分类条：每个子面板独立导出 ====
common<-theme_bw(base_family="Arial",base_size=8)+
  theme(panel.grid=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.spacing.x=unit(1.4,"mm"),strip.text=element_text(size=7),
        strip.background=element_rect(fill="#dedede",linewidth=.25),
        axis.text.y=element_text(colour="black",size=7),axis.title.y=element_text(size=7.5),
        panel.border=element_rect(colour="#666666",linewidth=.25),
        legend.title=element_text(size=7),legend.text=element_text(size=6.5),
        plot.margin=margin(1,1,1,1,"mm"))
p1<-ggplot(d,aes(TFcombo,CellState))+
  geom_tile(aes(fill=PercentActive))+
  geom_point(data=d[d$P_active<.05,],colour="white",size=.22)+
  scale_fill_gradientn(colours=c("black","red","orange"),limits=c(0,100),name="% sequences\nactive")+
  facet_grid(~Faceting,scales="free_x",space="free_x")+
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+labs(x=NULL,y="Cell state")+common+
  guides(fill=guide_colourbar(barheight=unit(17,"mm"),barwidth=unit(2.7,"mm")))
p2<-ggplot(d,aes(TFcombo,CellState))+
  geom_tile(aes(fill=PercentRepressed))+
  geom_point(data=d[d$P_repressed<.05,],colour="white",size=.22)+
  scale_fill_gradientn(colours=c("black","blue","green"),limits=c(0,35),oob=scales::squish,name="% sequences\nrepressed")+
  facet_grid(~Faceting,scales="free_x",space="free_x")+
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+labs(x=NULL,y="Cell state")+common+
  guides(fill=guide_colourbar(barheight=unit(17,"mm"),barwidth=unit(2.7,"mm")))
p3<-ggplot(summary,aes(TFcombo,1,fill=Class))+
  geom_tile()+facet_grid(~Faceting,scales="free_x",space="free_x")+
  scale_fill_manual(values=c(Activator="#66C2A5",Dual="#FC8D62",Repressor="#E78AC3",Other="#8DA0CB"),drop=FALSE,name="Factor-level\nsummary")+
  scale_x_discrete(expand=c(0,0))+scale_y_continuous(expand=c(0,0),breaks=NULL)+labs(x=NULL,y=NULL)+common+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=5.2,colour="black"),strip.text=element_blank(),strip.background=element_blank(),
        legend.key.size=unit(2.2,"mm"),legend.spacing.y=unit(.3,"mm"))
save_gg(p1,file.path(draft,"figure/08a1_active_panel"),180,42)
save_gg(p2,file.path(draft,"figure/08a2_repressed_panel"),180,42)
save_gg(p3,file.path(draft,"figure/08a3_class_panel"),180,43)

# 完整参考结构也保留一份；共享 gtable 列宽保证 70 个因子跨三行严格对齐。
g1<-ggplotGrob(p1);g2<-ggplotGrob(p2);g3<-ggplotGrob(p3)
widths<-unit.pmax(g1$widths,g2$widths,g3$widths)
g1$widths<-widths;g2$widths<-widths;g3$widths<-widths
g3$heights[unique(g3$layout$t[grepl("^panel",g3$layout$name)])]<-unit(2.2,"mm")
save_draw({
  grid.newpage();pushViewport(viewport(layout=grid.layout(3,1,heights=unit(c(36,36,25),"mm"))))
  pushViewport(viewport(layout.pos.row=1));grid.draw(g1);popViewport()
  pushViewport(viewport(layout.pos.row=2));grid.draw(g2);popViewport()
  pushViewport(viewport(layout.pos.row=3));grid.draw(g3);popViewport()
  popViewport()
},file.path(draft,"figure/08a_factor_state_heatmaps"),180,97)

## 3. 转录因子配对：上左红色激活、下右蓝绿抑制，白点标记 ====
pairs<-read.csv(file.path(module_root,"data/paired_factor_counts.csv"))
tfs<-unique(pairs$TF1)
pairs$TF1<-factor(pairs$TF1,levels=tfs);pairs$TF2<-factor(pairs$TF2,levels=tfs)
pairs$i<-as.integer(pairs$TF1);pairs$j<-as.integer(pairs$TF2)
pairs$PercentActive<-pairs$Active/pairs$N*100;pairs$PercentRepressed<-pairs$Repressed/pairs$N*100
pairs$P_active<-p.adjust(pbinom(pairs$Active-1,pairs$N,background_active,lower.tail=FALSE),"BH")
pairs$P_repressed<-p.adjust(pbinom(pairs$Repressed-1,pairs$N,background_repressed,lower.tail=FALSE),"BH")
write.csv(pairs,file.path(draft,"result/08b_paired_statistics.csv"),row.names=FALSE)
p4<-ggplot()+
  geom_tile(data=pairs[pairs$i<pairs$j,],aes(TF1,TF2,fill=PercentActive),colour="#353535",linewidth=.12)+
  scale_fill_gradientn(colours=c("black","red","orange","yellow"),limits=c(0,100),name="% active",
                       guide=guide_colourbar(barheight=unit(28,"mm"),barwidth=unit(3,"mm")))+
  new_scale_fill()+
  geom_tile(data=pairs[pairs$i>pairs$j,],aes(TF1,TF2,fill=PercentRepressed),colour="#353535",linewidth=.12)+
  scale_fill_gradientn(colours=c("black","blue","green"),limits=c(0,30),oob=scales::squish,name="% repressed",
                       guide=guide_colourbar(barheight=unit(28,"mm"),barwidth=unit(3,"mm")))+
  geom_tile(data=pairs[pairs$i==pairs$j,],aes(TF1,TF2),fill="white",colour="#777777",linewidth=.15)+
  geom_point(data=pairs[pairs$i<pairs$j & pairs$P_active<.05,],aes(TF1,TF2),colour="white",size=.25)+
  geom_point(data=pairs[pairs$i>pairs$j & pairs$P_repressed<.05,],aes(TF1,TF2),colour="white",size=.25)+
  scale_x_discrete(limits=tfs,position="top",expand=c(0,0))+scale_y_discrete(limits=tfs,expand=c(0,0))+
  coord_fixed()+labs(x=NULL,y=NULL)+theme_bw(base_family="Arial",base_size=8)+
  theme(panel.grid=element_blank(),axis.text.x=element_text(angle=90,hjust=0,vjust=.5,size=7,colour="black"),
        axis.text.y=element_text(size=7,colour="black"),axis.ticks=element_blank(),
        legend.title=element_text(size=7.5),legend.text=element_text(size=7),plot.margin=margin(3,3,3,3,"mm"))
save_gg(p4,file.path(draft,"figure/08b_paired_activity_heatmap"),180,180)

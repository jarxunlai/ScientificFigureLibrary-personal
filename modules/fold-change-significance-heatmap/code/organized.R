# 显著性星号倍数变化热图 / Fold-change heatmap with significance stars
# 科学问题：不同处理对代谢物的效应方向、幅度与显著性如何共同呈现？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "fold-change-significance-heatmap")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复纵向窄热图比例、类别分隔和星号可读性，保持FC色阶与量纲；输入/派生值均为合成示例，不冒充来源结果。
# 显著性星号热图；输入 log10FC 已经转换，不能再次取对数。
# 示例 P 来自生成器的示例 log10FC/SE；不是论文检验。星号阈值明确为 P<0.05。
library(ggplot2)
library(grid)


source(file.path(module_root,"code/export.R"))

## 1. 输入、顺序、类别分隔与星号 ====
d<-read.csv(file.path(module_root,"data/fold_changes.csv"))
d$BileAcid<-factor(d$BileAcid,levels=rev(unique(d$BileAcid)))
# 以原拟态数值恢复来源八个核受体标签；这不是作者实测响应。
receptors<-c("FXR","VDR","CAR-3","PXR","AHR","PPARA","PPARD","PPARG")
d$Treatment<-factor(receptors[match(d$Treatment,unique(d$Treatment))],levels=receptors)
d$Category<-factor(d$Category,levels=c("Conventional bile acids","BBAAs"))
d$Star<-ifelse(d$P<.05,"*","")
write.csv(d,file.path(draft,"result/05_plot_input.csv"),row.names=FALSE)

## 2. 按类别自由高度分面，保留青-白-珊瑚红与非均匀色阶 ====
p<-ggplot(d,aes(Treatment,BileAcid))+
  geom_tile(aes(fill=log10FC),colour="#555555",linewidth=.18)+
  geom_text(aes(label=Star,colour=log10FC<0),size=2.8,vjust=.7,show.legend=FALSE)+
  facet_grid(Category~.,scales="free_y",space="free_y",switch="y")+
  scale_fill_gradientn(colours=c("#6FCFCF","white","#DD6048"),
                       values=scales::rescale(c(-.35,0,2)),limits=c(-.35,2),
                       breaks=c(0,1,2),labels=c(expression(10^0),expression(10^1),expression(10^2)),name="Fold\nchange")+
  scale_colour_manual(values=c("FALSE"="#7C2115","TRUE"="#145444"))+
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))+
  labs(x=NULL,y=NULL)+
  theme_bw(base_family="Arial",base_size=8)+
  theme(panel.grid=element_blank(),panel.border=element_rect(fill=NA,colour="black",linewidth=.35),
        strip.placement="outside",strip.background=element_blank(),strip.text.y.left=element_text(angle=90,size=7.5),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=7),axis.text.y=element_text(colour="black",size=7.3),
        axis.ticks=element_line(linewidth=.2),axis.ticks.length=unit(1,"mm"),
        panel.spacing=unit(2,"mm"),legend.position="top",legend.title=element_text(size=7.5),legend.text=element_text(size=7),
        plot.margin=margin(2,3,2,2,"mm"))+
  guides(fill=guide_colourbar(barwidth=unit(42,"mm"),barheight=unit(3,"mm"),title.position="left",frame.colour="black"))
save_gg(p,file.path(draft,"figure/05_significance_heatmap"),89,180)

# 邻域丰度差异蜂群图 / Neighbourhood differential-abundance beeswarms
# 科学问题：不同细胞类型中的邻域丰度效应如何沿空间或时间对比变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "neighbourhood-abundance-beeswarm")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复邻域蜂群的可见点径、细胞类型行距及四对比布局，不运行MiloR上游；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 细胞丰度差异蜂群图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1011)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

library(ggbeeswarm)
df <- read.csv(file.path(module_root,"data/neighbourhoods.csv"))
# age-skeleton示例收敛到来源±1范围；其余比较维持独立尺度。
df$logFC[df$Comparison=="age_skeleton"]<-df$logFC[df$Comparison=="age_skeleton"]*.14
write.csv(df,file.path(base,"result/11_visual_neighbourhood_input.csv"),row.names=FALSE)
cell_order <- unique(df$cell_type)
df$cell_type <- factor(df$cell_type,levels=rev(cell_order))
panels <- list()
cols <- c(anterior_posterior="#f33131",age_skeleton="#4658a5",cranial_axis="#864093",age_cranial="#4658a5")
## 3. 使用预计算的邻域logFC绘制蜂群，不运行MiloR ====
for (contrast in names(cols)) {
  dd <- subset(df,Comparison==contrast)
  ref <- if(contrast=="age_skeleton")c(-1,-.5,.5,1) else c(-5,5)
  p <- ggplot(dd,aes(logFC,cell_type)) + geom_vline(xintercept=0,linetype=2,linewidth=.6) + geom_vline(xintercept=ref,colour="grey65",linewidth=.25) +
    geom_quasirandom(colour=cols[contrast],size=.48,width=.35,orientation="y",alpha=.7) +
    labs(x="logFC",y=NULL,title=c(anterior_posterior="Appendicular  <-->  Cranium",age_skeleton="5 PCW  <-->  11 PCW",cranial_axis="Anterior  <-->  Posterior",age_cranial="5 PCW  <-->  11 PCW")[contrast]) + theme_bw(base_size=7) + theme(panel.grid=element_blank(),axis.text.y=element_text(size=7))
  panels[[contrast]] <- p
  save_plot(p,file.path(fig_dir,paste0("11_",contrast)),88,142)
}
save_plot(wrap_plots(panels,ncol=2)+plot_annotation(tag_levels="a"),file.path(fig_dir,"11a_skeleton_beeswarm"),190,245)
## 4. 早期至晚期妊娠，连续色编码与谱系色条 ====
skin <- read.csv(file.path(module_root,"data/skin_neighbourhoods.csv"))
skin$cell_type <- factor(skin$cell_type,levels=rev(unique(skin$cell_type)))
stripe <- skin %>% distinct(cell_type,lineage_colour)
p <- ggplot(skin,aes(logFC,cell_type)) + geom_vline(xintercept=0,linetype=2,linewidth=.3,colour="grey55") +
  geom_quasirandom(aes(colour=logFC),size=.46,width=.35,orientation="y") + scale_colour_gradient2(low="#b47880",mid="white",high="#6659a3",midpoint=0,guide="none") +
  geom_point(data=stripe,aes(x=-1.2,y=cell_type,fill=lineage_colour),inherit.aes=FALSE,shape=22,size=2.2,colour="transparent") + scale_fill_identity() +
  coord_cartesian(xlim=c(-1.25,1.1),clip="off") + scale_x_continuous(breaks=c(-1,0,1)) +
  labs(x="Early gestation   <--- log(fold change) --->   Late gestation",y=NULL) + theme(axis.text.y=element_text(size=6))
save_plot(p,file.path(fig_dir,"11b_skin_beeswarm"),120,150)

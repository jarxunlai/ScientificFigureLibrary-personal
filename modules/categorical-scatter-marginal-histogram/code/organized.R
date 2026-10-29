# 分类散点、异方差误差线与边际直方图 / Categorical scatter with asymmetric errors and marginal histogram
# 科学问题：对象的两个测量及其不确定性是否随类别和时间发生变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "categorical-scatter-marginal-histogram")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：修正分类点形状为圆/方/菱形并恢复对数速度轴，恢复双向误差线和边缘直方图比例；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现| 带双线误差线的分类散点图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1004)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/runaways.csv"))
spectra <- c("WN(h)",paste0("O",2:9),"B or later","Unknown")
cols <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026","#808080")
df$spectral <- factor(df$spectral,levels=spectra)
df$spectral_id <- as.numeric(df$spectral)-1
## 3. 上方直方图/KDE，下方异方差双向误差线 ====
de <- density(df$t_kin,bw=.2,from=0,to=3)
kde <- data.frame(x=de$x,y=de$y)
hist_scale <- .15*nrow(df)
p_top <- ggplot(df,aes(t_kin)) + geom_histogram(binwidth=.15,boundary=0,fill="#f5a5a5",colour="white",linewidth=.2) +
  geom_line(data=kde,aes(x,y*hist_scale),inherit.aes=FALSE,colour="grey25",linewidth=.6) +
  scale_x_continuous(limits=c(0,3),expand=c(0,0)) + scale_y_continuous(sec.axis=sec_axis(~./hist_scale,name="Density")) + labs(y="Runaways",x=NULL) + theme_bw(base_size=8) + theme(panel.grid=element_blank())
p <- ggplot(df,aes(t_kin,v_t)) +
  geom_errorbar(aes(ymin=v_t-v_err,ymax=v_t+v_err),width=0,colour="grey55",linewidth=.35) +
  geom_errorbar(aes(xmin=t_kin-t_min,xmax=t_kin+t_plus),orientation="y",width=0,colour="grey55",linewidth=.35) +
  geom_point(aes(fill=spectral_id,shape=magnitude),size=1.9,stroke=.4) +
  scale_shape_manual(values=c(faint=21,intermediate=22,bright=23)) + scale_fill_gradientn(colours=cols,breaks=0:10,labels=spectra,limits=c(0,10),name="Spectral type",guide=guide_colourbar(barheight=unit(65,"mm"),barwidth=unit(3,"mm"),reverse=TRUE)) +
  scale_x_continuous(limits=c(0,3),expand=c(0,0)) + scale_y_log10(breaks=c(20,30,50,100,150,200,250),labels=scales::label_number(),expand=c(0,0)) + coord_cartesian(ylim=c(20,250)) +
  labs(x="Kinematic age (Myr)",y="Tangential velocity (km/s)") + theme_bw(base_size=8) + theme(panel.grid=element_blank(),legend.text=element_text(size=6))
save_plot(p_top,file.path(fig_dir,"04a_kinematic_age_histogram"),125,37)
save_plot(p,file.path(fig_dir,"04b_runaway_errors"),125,112)
save_plot(p_top / p + plot_layout(heights=c(.8,3.5)),file.path(fig_dir,"04_runaway_composite"),130,155)

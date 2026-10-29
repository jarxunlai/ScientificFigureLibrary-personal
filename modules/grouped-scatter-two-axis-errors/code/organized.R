# 分组散点、二次趋势与双向误差线 / Grouped scatter with two-axis error bars
# 科学问题：两个连续测量在不同组中如何共同变化，组内离散程度多大？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "grouped-scatter-two-axis-errors")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复来源寿命窗口及theory/五组/mean统一图例；输入/派生值均为合成示例，不冒充来源结果。
# 如何把普通的散点图，做成顶刊里的样子？
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1001)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/points.csv"))
cols <- c("7nm"="#4C9DCA","9nm"="#89CAF0","13nm"="#4EBB9E","16nm"="#EFBB4D","23nm"="#E28F4D")
df$Group <- factor(df$Group, levels = names(cols))
summary <- df %>% group_by(Group) %>% summarise(X=mean(Energy),Y=mean(Lifetime),X_sd=sd(Energy),Y_sd=sd(Lifetime))
write.csv(summary,file.path(root,"result/group_summary.csv"),row.names=FALSE)
## 3. 散点、二次趋势与均值双向标准差 ====
p <- ggplot(df,aes(Energy,Lifetime)) +
  geom_point(aes(colour=Group),size=1.15,alpha=.8) +
  geom_smooth(aes(linetype="theory"),method="lm",formula=y~poly(x,2),se=FALSE,colour="grey30",linewidth=.5) +
  geom_errorbar(data=summary,aes(x=X,ymin=Y-Y_sd,ymax=Y+Y_sd),inherit.aes=FALSE,width=.0015,linewidth=.4) +
  geom_errorbar(data=summary,aes(y=Y,xmin=X-X_sd,xmax=X+X_sd),inherit.aes=FALSE,orientation="y",width=12,linewidth=.4) +
  geom_point(data=summary,aes(X,Y,shape="mean"),inherit.aes=FALSE,size=2.5) +
  scale_colour_manual(values=cols) + scale_shape_manual(values=18,guide="none") + scale_linetype_manual(values="dashed",guide="none") +
  scale_x_continuous(breaks=seq(2.32,2.44,.02),sec.axis=dup_axis(name=NULL,labels=NULL)) +
  scale_y_continuous(breaks=seq(100,800,100),sec.axis=dup_axis(name=NULL,labels=NULL)) +
  labs(x="Energy (eV)",y="Lifetime (ps)") + coord_cartesian(xlim=c(2.31,2.45),ylim=c(50,850)) +
  theme_bw(base_size=9,base_family="Arial") + theme(panel.grid=element_blank(),legend.position="inside",legend.position.inside=c(.03,.97),legend.justification=c(0,1),legend.title=element_blank(),legend.spacing.y=unit(.3,"mm"),legend.key.height=unit(3,"mm"))
# 一份统一图例保留理论曲线及均值菱形，避免拆成三个重复的自动图例。
p<-p+theme(legend.position="none")+
 annotate("segment",x=2.316,xend=2.325,y=790,yend=790,linetype=2,linewidth=.45,colour="grey30")+
 annotate("text",x=2.329,y=790,label="theory",hjust=0,size=2.3)+
 annotate("point",x=rep(2.3205,5),y=seq(750,590,-40),colour=unname(cols),size=1.8)+
 annotate("text",x=rep(2.329,5),y=seq(750,590,-40),label=names(cols),hjust=0,size=2.3)+
 annotate("point",x=2.3205,y=550,shape=18,size=2.2)+
 annotate("text",x=2.329,y=550,label="mean",hjust=0,size=2.3)
save_plot(p,file.path(fig_dir,"01_energy_lifetime"),115,105)

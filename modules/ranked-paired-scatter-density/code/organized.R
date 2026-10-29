# 排序配对散点与边际密度 / Ranked paired scatter with marginal density
# 科学问题：同一对象在两个条件下的数值如何随差值排序而变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "ranked-paired-scatter-density")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复来源的细灰配对连接线，保持配对排序并对齐主图和右侧密度范围；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 两个条件下的配对散点图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1006)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/paired.csv")) %>% arrange(AB-no_AB) %>% mutate(Index=row_number())
long <- pivot_longer(df,c(AB,no_AB),names_to="Condition",values_to="Value")
cols <- c(AB="#223f9a",no_AB="#aaaaaa")
## 3. 每条竖线连接同一个triplet在两个条件下的值 ====
p <- ggplot() +
  geom_segment(data=df,aes(x=Index,xend=Index,y=AB,yend=no_AB),colour="grey80",linewidth=.14) +
  geom_smooth(data=long,aes(Index,Value,colour=Condition,fill=Condition),method="loess",formula=y~x,linewidth=.55,se=TRUE) +
  geom_point(data=long,aes(Index,Value,fill=Condition),shape=21,colour="grey40",size=1.0,stroke=.15) +
  scale_colour_manual(values=cols,labels=c("BC given AB coordination","BC given no AB coordination")) + scale_fill_manual(values=cols,labels=c("BC given AB coordination","BC given no AB coordination")) +
  coord_cartesian(ylim=c(-2.3,3.5)) + labs(x="all triplets (sorted on increasing delta Z)",y="BC coordination") +
  annotate("segment",x=110,xend=135,y=-2,yend=-2,arrow=arrow(type="closed",length=unit(1.3,"mm")),linewidth=.3) + annotate("text",x=140,y=-2,label="stronger triplets",hjust=0,size=2.4) +
  theme(legend.position="inside",legend.position.inside=c(.02,.99),legend.justification=c(0,1),legend.text=element_text(size=7))
pd <- ggplot(long,aes(y=Value,fill=Condition)) + geom_density(aes(x=after_stat(density)),orientation="y",alpha=.5,linewidth=.25) + scale_fill_manual(values=cols) +
  coord_cartesian(ylim=c(-2.3,3.5),expand=FALSE) + theme_void() + theme(legend.position="none",plot.margin=margin(5,0,5,0))
save_plot(p,file.path(fig_dir,"06a_paired_sorted"),125,95)
save_plot(pd,file.path(fig_dir,"06b_marginal_density"),25,95)
save_plot(p+pd+plot_layout(widths=c(5,.55)),file.path(fig_dir,"06_paired_with_density"),140,95)

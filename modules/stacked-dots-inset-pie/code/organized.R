# 堆积病灶点图与内嵌饼图 / Stacked lesion dots with an inset summary pie
# 科学问题：每个对象的观测是否被采样覆盖，整体覆盖比例是多少？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "stacked-dots-inset-pie")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复堆积病灶的圆点大小与病例标签，内嵌饼图仍由相同病灶表聚合；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 堆积点图 + 内嵌饼图
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1022)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/metastases.csv"))
cols <- c("imaged + sampled"="#c8102e","sampled only"="#e6b794","imaged only"="white")
counts <- df %>% count(id) %>% arrange(n)
df$id<-factor(df$id,levels=counts$id);df$group<-factor(df$group,levels=names(cols))
df <- df %>% arrange(id,group) %>% group_by(id) %>% mutate(group_seq=row_number()) %>% ungroup()
## 3. 每个点一处转移灶，数量排序与原图结构一致 ====
pm <- ggplot(df,aes(id,group_seq,fill=group)) + geom_point(shape=21,colour="black",stroke=.3,size=1.8) + scale_fill_manual(values=cols) +
  scale_y_continuous(limits=c(0,32),expand=c(0,0),breaks=seq(0,30,5)) + labs(x=NULL,y="Metastases") + theme(axis.text.x=element_text(angle=90,size=6.8,hjust=1,vjust=.5),legend.position="inside",legend.position.inside=c(.02,.99),legend.justification=c(0,1),legend.text=element_text(size=6))
pie <- df %>% count(group)
pp <- ggplot(pie,aes(x=1,y=n,fill=group)) + geom_col(width=1,colour="black",linewidth=.4) + coord_polar(theta="y",start=0) +
  geom_text(aes(label=n),position=position_stack(vjust=.5),size=3) + scale_fill_manual(values=cols) + theme_void() + theme(legend.position="none",plot.margin=margin(2,2,2,2))
save_plot(pm,file.path(fig_dir,"22a_stacked_dots"),150,100);save_plot(pp,file.path(fig_dir,"22b_summary_pie"),65,65)
save_plot(pm+inset_element(pp,left=.18,bottom=.46,right=.58,top=.94),file.path(fig_dir,"22_stacked_dot_inset_pie"),150,105)
write.csv(pie,file.path(root,"result/summary_counts.csv"),row.names=FALSE)

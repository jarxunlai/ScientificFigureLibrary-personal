# 气泡回归与肿瘤富集关联散点图 / Bubble regressions and tumor-enrichment association scatter
# 科学问题：两种测量是否一致，富集程度与疗效指标有何关联？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "bubble-regression-enrichment-scatter")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：纠正PIRATA站点颜色、气泡面积编码；补齐正负两种肿瘤反应图，不混淆方向；输入/派生值均为合成示例，不冒充来源结果。
# 顶刊复现 | 气泡图/散点图+线性回归
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1007)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/ocean.csv"))
cols <- c("KEO"="#07ddfb","RAMAS"="#fb5454","WHOTS"="#179dfc","IMOS"="#fc7c1f","MOVE3"="#8484fc","Stratus"="#fcbc10","RAMAN"="#2e3afc","PIRATA"="#fcfc00")
panels <- list()
## 3. 四个观测比较面板；后二者为过原点的回归 ====
for (metric in c("EHTA","ELTA","pEHTA","pELTA")) {
  dd <- subset(df,Metric==metric)
  p <- ggplot(dd,aes(x,y)) + geom_abline(slope=1,intercept=0,linetype=2,linewidth=.4) + geom_point(aes(colour=Site,size=Depth),alpha=.8) +
    scale_colour_manual(values=cols) + scale_size_area(max_size=5.0,breaks=c(10,100,1000),labels=c("10 m","100 m","1000 m")) + labs(x=paste(metric,"mooring"),y=paste(metric,"WOD")) + coord_fixed() + theme(legend.position="bottom",legend.text=element_text(size=6))
  if (metric %in% c("pEHTA","pELTA")) {
    fit <- lm(y~0+x,data=dd); sm <- summary(fit)
    p <- p + geom_smooth(method="lm",formula=y~x+0,colour="black",fill="grey80",linewidth=.5) +
      annotate("text",x=0,y=80,label=sprintf("y = %.2f x; adj. R2 = %.2f\nP = %.2g",coef(fit)[1],sm$adj.r.squared,sm$coefficients[1,4]),hjust=0,vjust=1,size=2.1)
  }
  panels[[metric]] <- p
  save_plot(p,file.path(fig_dir,paste0("07_",metric)),105,105)
}
save_plot((wrap_plots(panels,ncol=2,guides="collect")+plot_annotation(tag_levels="a")) & theme(legend.position="bottom"),file.path(fig_dir,"07a_ocean_four_panel"),180,170)
## 4. 肿瘤富集与疗效关联；保留负相关的符号 ====
dd <- read.csv(file.path(module_root,"data/tumor_response.csv"))
dd$Ti_original<-dd$Ti;dd$Ti<--dd$Ti
fit <- lm(Ti~Tumor,data=dd); rr <- cor.test(dd$Tumor,dd$Ti)
p <- ggplot(dd,aes(Tumor,Ti)) + geom_smooth(method="lm",formula=y~x,colour="black",fill="grey80",linewidth=.5) +
  geom_point(aes(colour=SubCellType),size=1.8) + geom_hline(yintercept=0,linetype=3,linewidth=.3) + geom_vline(xintercept=0,linetype=3,linewidth=.3) +
  geom_text_repel(data=subset(dd,label==1),aes(label=SubCellType),size=2,min.segment.length=0,max.overlaps=Inf,seed=47) +
  annotate("text",x=-90,y=.5,label=sprintf("r = %.2f\nP = %.2g",unname(rr$estimate),rr$p.value),hjust=0,vjust=1,size=2.5) +
  scale_colour_manual(values=scales::hue_pal()(nrow(dd))) + labs(x="Tumor enrichment\nEnriched in normal  <---->  Enriched in tumor",y="Better response (Ti)\nGreater decrease  <---->  Greater increase") + theme(legend.position="none")
save_plot(p,file.path(fig_dir,"07b_tumor_response"),130,115)

# 原文另有正向疗效版本；作为独立完整图保存，而不是替换负向语义。
dd$Ti<-dd$Ti_original
p_positive<-ggplot(dd,aes(Tumor,Ti))+geom_smooth(method="lm",colour="black",fill="grey85",linewidth=.5)+
 geom_point(aes(colour=SubCellType),size=2)+geom_text_repel(data=subset(dd,label==1),aes(label=SubCellType),size=2.5,seed=47,max.overlaps=Inf)+
 labs(x="Tumor enrichment",y="Response increase (illustrative)")+theme(legend.position="none")
save_plot(p_positive,file.path(fig_dir,"07c_tumor_response_increase"),130,115)
write.csv(dd,file.path(base,"result/07_tumor_visual_input.csv"),row.names=FALSE)

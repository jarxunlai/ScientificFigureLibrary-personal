# 十种相关性图复刻案例总览 / Correlation case-study collection
# 科学问题：不同变量维度与数据结构应选用何种相关性表达？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "correlation-figure-case-collection")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复对数刻度、MHQ类别图例、图内残差框及多面板标识；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure/case-corr06")
dat <- file.path(module_root,"data/case-corr06")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 基因组大小与编码序列数；残差在log尺度计算 ====
d <- read.csv(file.path(dat,"genomes.csv"))
d$type <- factor(d$type,levels=c("SAG","MAG","WGS"))
fit <- lm(log(number_of_cds)~log(length),data=d); d$residual <- residuals(fit)
write.csv(d,file.path(root,"result/case-corr06","genomes_with_residuals.csv"),row.names=FALSE)
cols <- c(SAG="#1F77B4",MAG="#FF7F0E",WGS="#808080")
## 2. 主散点、上下同源边缘直方图、图内残差箱线图 ====
p <- ggplot(d,aes(length,number_of_cds,colour=type))+
  geom_point(size=.95,alpha=.65)+geom_abline(intercept=coef(fit)[1]/log(10),slope=coef(fit)[2],linetype=2,colour="#ef5a4d",linewidth=.5)+
  scale_x_log10(breaks=c(1e6,1e7),labels=scales::label_log())+scale_y_log10(breaks=c(1e3,1e4),labels=scales::label_log())+scale_colour_manual(values=cols,labels=paste0("MHQ ",names(cols)," (N=",as.integer(table(d$type)),")"),name=NULL)+
  labs(x="Genome total size [bp] (corrected)",y="Number of predicted CDS (corrected)")+
  theme(legend.position="inside",legend.position.inside=c(.03,.98),legend.justification=c(0,1),legend.text=element_text(size=7),panel.border=element_rect(fill=NA,colour="grey60",linewidth=.35))
top <- ggplot(d,aes(log10(length),fill=type))+geom_histogram(bins=55,position="identity",alpha=.28,colour="grey50",linewidth=.15)+scale_fill_manual(values=cols)+theme_void()+theme(legend.position="none")
right <- ggplot(d,aes(log10(number_of_cds),fill=type))+geom_histogram(bins=55,position="identity",alpha=.28,colour="grey50",linewidth=.15)+scale_fill_manual(values=cols)+coord_flip()+theme_void()+theme(legend.position="none")
box <- ggplot(d,aes(type,residual,fill=type))+geom_boxplot(width=.6,outlier.size=.35,linewidth=.25)+scale_fill_manual(values=cols)+coord_flip()+labs(title="Residuals (log)",x=NULL,y=NULL)+theme_bw(base_size=7)+theme(legend.position="none",panel.grid=element_blank(),plot.title=element_text(size=8,hjust=.5),plot.margin=margin(1,1,1,1))
p <- p+inset_element(box,left=.49,bottom=.04,right=.96,top=.28)
out <- (top+plot_spacer()+plot_layout(widths=c(1,.19)))/(p+right+plot_layout(widths=c(1,.19)))+plot_layout(heights=c(.19,1))
out<-out+plot_annotation(tag_levels=list(c("b")))
save_plot(out,file.path(fig,"06_scatter_marginal_hist_residual"),175,170)

# 三种平滑散点曲线 / Three smooth scatter-curve patterns
# 科学问题：剂量响应、单峰变化与光谱峰形分别适用什么平滑方法？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "smooth-scatter-curve-patterns")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：恢复浓度对数刻度和两类点形，强化电导峰值与平滑光谱的几何层级；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 散点图+平滑曲线（3种方法）
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1005)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

dose <- read.csv(file.path(module_root,"data/dose_response.csv"))
conductivity <- read.csv(file.path(module_root,"data/conductivity.csv"))
spectra <- read.csv(file.path(module_root,"data/spectra.csv"))
## 3. 四参数logistic拟合：参数与输入同文件保留 ====
dose_fits <- list()
for (g in unique(dose$condition)) {
  dd <- subset(dose,condition==g)
  fit <- nls(shift~SSfpl(log10(concentration),A,B,C,D),data=dd)
  xx <- 10^seq(-9,-4,length.out=300)
  dose_fits[[g]] <- data.frame(concentration=xx,shift=predict(fit,newdata=data.frame(concentration=xx)),condition=g)
}
fit_df <- bind_rows(dose_fits)
p1 <- ggplot(dose,aes(concentration,shift,colour=condition)) + geom_line(data=fit_df,linewidth=.5) + geom_point(aes(shape=condition,fill=condition),colour="grey25",size=2.0,stroke=.4) +
  scale_colour_manual(values=c("No Venetoclax"="#55a0fb",Venetoclax="#ff8000")) + scale_fill_manual(values=c("No Venetoclax"="#55a0fb",Venetoclax="#ff8000")) + scale_shape_manual(values=c(22,21)) +
  scale_x_log10(breaks=10^seq(-9,-4),labels=scales::label_log()) + labs(x="Concentration (M)",y="Shift (nm)",title="DBVen1619_2") + theme(legend.position="inside",legend.position.inside=c(.02,.98),legend.justification=c(0,1),legend.text=element_text(size=7))
save_plot(p1,file.path(fig_dir,"05a_logistic_dose_response"),110,95)
## 4. 形状保持样条；峰值渐变圈只作视觉强调 ====
fit2 <- list()
for (g in unique(conductivity$Group)) {
  dd <- subset(conductivity,Group==g); xx <- seq(min(dd$Concentration),max(dd$Concentration),length.out=300)
  fit2[[g]] <- data.frame(Concentration=xx,Conductivity=splinefun(dd$Concentration,dd$Conductivity,method="monoH.FC")(xx),Group=g)
}
peaks <- conductivity %>% group_by(Group) %>% slice_max(Conductivity,n=1)
cols <- c("LiFSI/FAN"="#2A328D","LiTFSI/FAN"="#0C8040","LiPF6/EC-EMC"="#E41E26")
p2 <- ggplot(conductivity,aes(Concentration,Conductivity,colour=Group))
for (s in seq(12,4,length.out=18)) p2 <- p2 + geom_point(data=peaks,size=s,alpha=.018,show.legend=FALSE)
p2 <- p2 + geom_line(data=bind_rows(fit2),linewidth=.55) + geom_point(shape=21,fill="white",size=1.9,stroke=.45) + scale_colour_manual(values=cols) +
  coord_cartesian(xlim=c(0,4.5),ylim=c(0,50)) + labs(x="Concentration (M)",y="Ionic conductivity (mS/cm)") + theme(legend.position="inside",legend.position.inside=c(.98,.98),legend.justification=c(1,1),legend.text=element_text(size=7))
save_plot(p2,file.path(fig_dir,"05b_conductivity_spline"),110,95)
## 5. 平滑样条拟合光谱；拟合不解释为实验结果 ====
fit3 <- list()
for (g in unique(spectra$group)) {
  dd <- subset(spectra,group==g); fit <- smooth.spline(dd$wavelength,dd$intensity,spar=.5)
  xx <- seq(min(dd$wavelength),max(dd$wavelength),length.out=1000)
  fit3[[g]] <- data.frame(wavelength=xx,intensity=predict(fit,xx)$y,group=g)
}
p3 <- ggplot(spectra,aes(wavelength,intensity,colour=group)) + geom_point(size=.25,alpha=.55) + geom_line(data=bind_rows(fit3),linewidth=.65) +
  scale_colour_manual(values=c("400 nN"="#ff7b0b","0 nN"="#27a9e0")) + labs(x="Wavelength (nm)",y="Normalized emission intensity") + theme(legend.position="inside",legend.position.inside=c(.05,.97),legend.justification=c(0,1))
save_plot(p3,file.path(fig_dir,"05c_spectral_smoothing"),100,100)
write.csv(fit_df,file.path(root,"result/dose_fitted.csv"),row.names=FALSE)
write.csv(bind_rows(fit2),file.path(root,"result/conductivity_fitted.csv"),row.names=FALSE)
write.csv(bind_rows(fit3),file.path(root,"result/spectra_fitted.csv"),row.names=FALSE)

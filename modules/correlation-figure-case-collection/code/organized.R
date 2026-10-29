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
# 子案例有独立结果目录；从空输出目录运行时也完整创建，不依赖本地旧产物。
for (case_id in c("case-corr01", "case-corr02", "case-corr03", "case-corr05", "case-corr06", "case-corr07")) {
  for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir, case_id), showWarnings = FALSE, recursive = TRUE)
}
# 第二轮：保留11个完整独立案例并复制其PDF配对，停止生成总拼图和preview替代图；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 十案例文章有11个图片块；完整保留，8项复用本批图，另补3类 ====
library(cowplot)
energy <- read.csv(file.path(dat,"energy_lifetime.csv"))
energy$Size <- factor(energy$Size,levels=c("7 nm","9 nm","13 nm","16 nm","23 nm"))
colours <- c("#4a9fc4","#8acde9","#45b69b","#f0bf42","#e48b47")
means <- energy %>% group_by(Size) %>% summarise(x=mean(Energy),y=mean(Lifetime),sx=sd(Energy),sy=sd(Lifetime),.groups="drop")
p2 <- ggplot(energy,aes(Energy,Lifetime,colour=Size))+geom_point(size=1.2,alpha=.85)+
  geom_abline(slope=4500,intercept=-10250,linetype=2,colour="grey30",linewidth=.4)+
  geom_errorbar(data=means,aes(x=x,y=y,ymin=y-sy,ymax=y+sy),inherit.aes=FALSE,width=.002,linewidth=.3)+
  geom_errorbar(data=means,aes(x=x,y=y,xmin=x-sx,xmax=x+sx),inherit.aes=FALSE,orientation="y",width=12,linewidth=.3)+
  geom_point(data=means,aes(x,y),inherit.aes=FALSE,shape=18,size=2.8)+
  scale_colour_manual(values=colours,name=NULL)+labs(x="Energy (eV)",y="Lifetime (ps)")+theme_bw(base_size=8)+
  theme(panel.grid=element_blank(),legend.position="inside",legend.position.inside=c(.02,.98),legend.justification=c(0,1),legend.key.height=unit(3,"mm"))
save_plot(p2,file.path(fig,"08_case02_energy_lifetime"),125,105)
## 2. 配对海洋指标：保留4个面板，每个也单独输出 ====
ocean <- read.csv(file.path(dat,"paired_ocean.csv"))
site_cols <- c(KEO="#00cee5",WHOTS="#159fff",MOVE3="#8d82ff",RAMAN="#333cff",RAMAS="#ff4d5b",IMOS="#ff7908",Stratus="#ffbd00",PIRATA="#e8e500")
panels <- list()
for(k in 1:4) {
  dd <- ocean[ocean$panel==k,]
  pp <- ggplot(dd,aes(x,y))+geom_point(aes(colour=Site,size=depth),alpha=.85)+
    scale_colour_manual(values=site_cols)+scale_size_area(max_size=3,breaks=c(10,100,1000),name="Depth (m)")+
    labs(x=if(k<=2) "Mooring anomaly (deg C)" else "Mooring proportion (%)",y=if(k<=2) "WOD anomaly (deg C)" else "WOD proportion (%)")+
    theme_classic(base_size=8)+theme(legend.position="none")
  if(k<=2) pp <- pp+geom_abline(slope=1,intercept=0,linetype=2,colour="grey35",linewidth=.4) else {
    fit <- lm(y~x,data=dd)
    pp <- pp+geom_smooth(method="lm",se=TRUE,colour="grey20",fill="grey75",linewidth=.5)+
      annotate("text",x=-Inf,y=Inf,label=sprintf("y = %.2fx %+.1f; R^2 = %.2f",coef(fit)[2],coef(fit)[1],summary(fit)$r.squared),hjust=-.04,vjust=1.2,size=2)
  }
  panels[[k]] <- pp
  save_plot(pp,file.path(fig,paste0("08_case03_ocean_panel",k)),115,105)
}
legend_plot <- panels[[1]]+theme(legend.position="bottom")+
  guides(colour=guide_legend(nrow=2,order=1,override.aes=list(size=2)),size=guide_legend(order=2))
legend_grob <- cowplot::get_legend(legend_plot)
ocean_grid <- cowplot::plot_grid(wrap_plots(panels,ncol=2),legend_grob,ncol=1,rel_heights=c(1,.14))
save_plot(ocean_grid,file.path(fig,"08_case03_paired_ocean"),190,175)
## 3. 配对协调图：每对点用浅线绑定；右侧同源密度 ====
co <- read.csv(file.path(dat,"coordination_triplets.csv"))
pc <- ggplot(co,aes(rank,with_AB))+
  geom_point(aes(colour="Given AB coordination"),size=.6)+geom_point(aes(y=no_AB,colour="No AB coordination"),shape=21,fill="grey90",size=.6,stroke=.25)+
  geom_smooth(aes(colour="Given AB coordination"),method="loess",se=TRUE,fill="#9ba8cf",linewidth=.5)+
  geom_smooth(aes(y=no_AB,colour="No AB coordination"),method="loess",se=TRUE,fill="grey80",linewidth=.5)+
  scale_colour_manual(values=c(`Given AB coordination`="#283e86",`No AB coordination`="grey60"),name=NULL)+
  labs(x="All triplets (sorted on increasing delta Z)",y="BC coordination")+theme_classic(base_size=8)+
  theme(legend.position="inside",legend.position.inside=c(.03,.98),legend.justification=c(0,1),legend.text=element_text(size=6))
long <- pivot_longer(co,cols=c(with_AB,no_AB),names_to="condition",values_to="coordination")
pd <- ggplot(long,aes(coordination,fill=condition))+geom_density(alpha=.5,linewidth=.3)+scale_fill_manual(values=c(with_AB="#8a9bc9",no_AB="grey80"))+coord_flip()+theme_void()+theme(legend.position="none")
save_plot(pc+pd+plot_layout(widths=c(1,.14)),file.path(fig,"08_case04_paired_coordination"),190,135)
## 4. 从本条目内的六份绘图实现实际重绘，不依赖其他Draft或预制截图 ====
sys.source(file.path(module_root,"code/case-corr01.R"),envir=new.env(parent=globalenv()))
sys.source(file.path(module_root,"code/case-corr02.R"),envir=new.env(parent=globalenv()))
sys.source(file.path(module_root,"code/case-corr03.R"),envir=new.env(parent=globalenv()))
sys.source(file.path(module_root,"code/case-corr05.R"),envir=new.env(parent=globalenv()))
sys.source(file.path(module_root,"code/case-corr06.R"),envir=new.env(parent=globalenv()))
sys.source(file.path(module_root,"code/case-corr07.R"),envir=new.env(parent=globalenv()))
paths <- c(
  file.path(batch,"figure/case-corr06/06_scatter_marginal_hist_residual.png"),
  file.path(fig,"08_case02_energy_lifetime.png"),
  file.path(fig,"08_case03_paired_ocean.png"),
  file.path(fig,"08_case04_paired_coordination.png"),
  file.path(batch,"figure/case-corr01/01_rice_matrix_year_trends.png"),
  file.path(batch,"figure/case-corr07/07_grouped_scatter_matrix_network.png"),
  file.path(batch,"figure/case-corr03/03_bile_acid_matrix_marginals.png"),
  file.path(batch,"figure/case-corr05/05_inverted_triangle_dendrogram.png"),
  file.path(batch,"figure/case-corr02/02a_grouped_environment_correlation.png"),
  file.path(batch,"figure/case-corr02/02b_cell_factor_clustered_heatmap.png"),
  file.path(batch,"figure/case-corr02/02c_mantel_environment_network.png")
)
# 原文章是11个独立案例，不再将它们拼成一张preview当作单图。
for(i in seq_along(paths)) {
 stem<-file.path(fig,sprintf("08_case%02d_complete",i))
 file.copy(paths[i],paste0(stem,".png"),overwrite=TRUE)
 file.copy(sub("\\.png$",".pdf",paths[i]),paste0(stem,".pdf"),overwrite=TRUE)
}
write.csv(data.frame(case=seq_along(paths),render=basename(paths)),file.path(entry,"result/08_case_mapping.csv"),row.names=FALSE)

# 分组散点与线性对数拟合六面板 / Grouped scatter with linear and logarithmic fits
# 科学问题：不同生态系统中初始存量与损失、干扰时间与相对损失有什么关联？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "grouped-scatter-log-linear-fit")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：修复拟合统计文字在各量纲中的放置，恢复三生态系统上下六面板与共同图例；输入/派生值均为合成示例，不冒充来源结果。
# Nature 复现 | 分组散点+线性/对数拟合线
# 目的：复建来源教程的绘图结构；输入为本条目 CSV 示例，不运行上游实验/组学分析。
# 输入：data/*.csv；输出：figure/*.png + *.pdf（300 dpi PNG / Cairo PDF）。
## 1. 环境与路径 ====
Sys.setlocale("LC_CTYPE", "English_United States.utf8")




source(file.path(module_root,"code/export.R"), encoding = "UTF-8")
library(ggrepel)
set.seed(1012)
fig_dir <- file.path(base, "figure")
theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_blank(), legend.key.height = unit(3.5, "mm"),
        plot.margin = margin(5, 6, 5, 5)))
## 2. 输入与可调参数 ====

df <- read.csv(file.path(module_root,"data/blue_carbon.csv"))
# 原CSV作为合成种子；按各生态系统恢复来源不同扰动分组/范围，不沿用三列完全一样的模型。
set.seed(271012)
df$Disturbance<-as.character(df$Disturbance)
source_groups<-list(Mangrove=c("Agriculture","Aquaculture","Climate/hydrological change","Harvesting"),
 Saltmarsh=c("Agriculture","Climate/hydrological change","Grazing"),Seagrass=c("Dredging","Vegetation cover damage"))
source_colors<-list(Mangrove=c("#ef9f19","#b51a22","#f329cd","#8e468f"),Saltmarsh=c("#278b77","#39c3c7","#7c891c"),Seagrass=c("#2922e6","#45aeea"))
visual<-list()
for(ecosystem in names(source_groups)){
 dd<-subset(df,Ecosystem==ecosystem);ng<-length(source_groups[[ecosystem]])
 dd$Disturbance<-source_groups[[ecosystem]][(seq_len(nrow(dd))-1)%%ng+1]
 dd$Initial<-dd$Initial/max(dd$Initial)*c(Mangrove=900,Saltmarsh=310,Seagrass=460)[ecosystem]
 dd$Duration<-dd$Duration/max(dd$Duration)*c(Mangrove=48,Saltmarsh=460,Seagrass=90)[ecosystem]+.3
 if(ecosystem=="Mangrove"){
   slope<-c(Agriculture=-.9,Aquaculture=-.78,"Climate/hydrological change"=-.55,Harvesting=-.2)
   dd$Absolute<-dd$Initial*slope[dd$Disturbance]+rnorm(nrow(dd),0,100)
   dd$Relative<-50-36*log(dd$Duration)+rnorm(nrow(dd),0,48)
 }else if(ecosystem=="Saltmarsh"){
   slope<-c(Agriculture=.6,"Climate/hydrological change"=-.6,Grazing=.25)
   dd$Absolute<-dd$Initial*slope[dd$Disturbance]+rnorm(nrow(dd),0,45)
   dd$Relative<-110-20*log(dd$Duration)+rnorm(nrow(dd),0,45)
 }else{
   dd$Absolute<-190-.85*dd$Initial+rnorm(nrow(dd),0,70)
   dd$Relative<-28-17*log(dd$Duration)+rnorm(nrow(dd),0,17)
 }
 visual[[ecosystem]]<-dd
}
df<-bind_rows(visual)
write.csv(df,file.path(base,"result/12_visual_synthetic_input.csv"),row.names=FALSE)
panels<-list()
## 3. 三种生态系统独立扰动、量纲和拟合；R平方/P均重新计算 ====
for(ecosystem in names(source_groups)){
 dd<-subset(df,Ecosystem==ecosystem)
 dd$Disturbance<-factor(dd$Disturbance,levels=source_groups[[ecosystem]])
 cols<-setNames(source_colors[[ecosystem]],source_groups[[ecosystem]])
 for(relation in c("absolute","relative")){
   if(relation=="absolute"){
     dd$x<-dd$Initial;dd$y<-dd$Absolute;model<-y~x
     xl<-"Initial soil C stock (Mg C ha-1)";yl<-"Soil C stock change (Mg C ha-1)"
   }else{
     dd$x<-dd$Duration;dd$y<-dd$Relative;model<-y~log(x)
     xl<-"Duration (Years)";yl<-"Soil C stock change (%)"
   }
   stats<-list()
   for(g in levels(dd$Disturbance)){
     ss<-subset(dd,Disturbance==g);fit<-lm(model,data=ss)
     stats[[g]]<-data.frame(Disturbance=g,label=sprintf("R2 = %.2f; P = %.2g",summary(fit)$r.squared,coef(summary(fit))[2,4]))
   }
   st<-bind_rows(stats);st$x<-min(dd$x);st$y<-max(dd$y)-seq(0,.20,length.out=nrow(st))*diff(range(dd$y))
   p<-ggplot(dd,aes(x,y,colour=Disturbance))+geom_hline(yintercept=0,linetype=3,linewidth=.3,colour="grey65")+
     geom_point(size=1.15,alpha=.8)+geom_smooth(method="lm",formula=model,linewidth=.5,se=TRUE,alpha=.1)+
     scale_colour_manual(values=cols)+geom_text(data=st,aes(x,y,label=label),inherit.aes=FALSE,colour=unname(cols),hjust=0,size=2.25)+
     labs(title=if(relation=="absolute")ecosystem else NULL,x=xl,y=yl)+theme_bw(base_size=8)+
     theme(panel.grid=element_blank(),legend.position=if(relation=="absolute")"inside"else"none",legend.position.inside=c(.98,.98),
       legend.justification=c(1,1),legend.text=element_text(size=5.5),legend.key.height=unit(2.4,"mm"),legend.key.width=unit(3,"mm"),
       legend.background=element_blank(),plot.title=element_text(size=9,hjust=.5),axis.title=element_text(size=7))
   if(relation=="absolute")p<-p+geom_abline(slope=-.25,intercept=0,linetype=2,colour="grey65",linewidth=.3)+geom_abline(slope=-.75,intercept=0,linetype=2,colour="grey35",linewidth=.3)
   panels[[paste(ecosystem,relation,sep="_")]]<-p
   save_plot(p,file.path(fig_dir,paste0("12_",ecosystem,"_",relation)),110,95)
 }
}
# 标签按原来源的a/b、c/d、e/f纵向配对，完整六面板不拆条目。
order<-c(1,3,5,2,4,6)
for(k in seq_along(order))panels[[order[k]]]<-panels[[order[k]]]+labs(tag=c("a","c","e","b","d","f")[k])
save_plot(wrap_plots(panels[order],ncol=3),file.path(fig_dir,"12_blue_carbon_six_panel"),220,155)

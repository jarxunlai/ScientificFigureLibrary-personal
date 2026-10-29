# 五类三元图与状态概率密度 / Ternary plots and state densities
# 科学问题：一个细胞、OTU或样本如何在三种相对状态或环境间分配？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "ternary-state-density")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：补回早期胚胎三元示例、恢复三元点大小层级，五类图仍独立输出；输入/派生值均为合成示例，不冒充来源结果。
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 三元比例保持非负、每行归一化；明确使用合成例子 ====
library(ggtern)
emb <- read.csv(file.path(dat,"embryo_lineages.csv"))
p1 <- ggtern(emb,aes(EPI,TE,PrE,colour=CellType))+geom_point(size=1.05,alpha=.75)+
  scale_colour_manual(values=c(EPI="#2273ae",PrE="#2c9e2c",TE="#fa7e0d"))+
  labs(title="36 h post iEFCs",x="EPI",y="TE",z="PrE")+theme_rgbg(base_size=9)+theme(legend.position="none",plot.title=element_text(hjust=.5))
save_plot(p1,file.path(fig,"10a_embryo_lineages"),150,150)
for(stage in c("BCl0_0h","iEFC0_0h")) {
 dd<-emb
 dd$CellType<-if(stage=="BCl0_0h")"EPI"else"PrE"
 dd$EPI<-if(stage=="BCl0_0h").75+.2*dd$EPI else .12+.15*dd$EPI
 dd$TE<-.05+.15*dd$TE;dd$PrE<-1-dd$EPI-dd$TE
 pp<-ggtern(dd,aes(EPI,TE,PrE,colour=CellType))+geom_point(size=.9,alpha=.6)+
  scale_colour_manual(values=c(EPI="#2273ae",PrE="#c561a6",TE="#fa7e0d"))+theme_rgbg(base_size=10)+
  labs(title=stage,x="EPI",y="TE",z="PrE")+theme(legend.position="none",plot.title=element_text(hjust=.5))
 save_plot(pp,file.path(fig,paste0("10a_",stage)),145,145)
 write.csv(dd,file.path(entry,"result",paste0("10a_",stage,"_visual.csv")),row.names=FALSE)
}
## 2. 免疫状态：年龄组离散色与年龄连续色分别出图 ====
im <- read.csv(file.path(dat,"immune_states.csv"))
for(mode in c("group","age")) {
  p <- ggtern(im,aes(TH2,TH1,TH17))+geom_point(aes(fill=if(mode=="group") AgeGroup else Age),shape=21,size=1.6,stroke=.15,colour="grey25",alpha=.85)+
    labs(x="TH2 state",y="TH1 state",z="TH17 state",title="CD4 TCM cell",fill=if(mode=="group") NULL else "Age (years)")+theme_rgbw(base_size=8)+
    theme(plot.title=element_text(hjust=.5),legend.position="right")
  if(mode=="group") p <- p+scale_fill_manual(values=c(Young="#35978f",Older="#bf812d")) else p <- p+scale_fill_distiller(palette="Spectral",direction=-1)
  save_plot(p,file.path(fig,paste0("10b_immune_state_",mode)),175,150)
}
## 3. 生态位富集：点大小为平均相对丰度 ====
otu <- read.csv(file.path(dat,"otu_niches.csv"))
p3 <- ggtern(otu,aes(Root,Nodule,Rhizosphere))+geom_point(aes(size=mean_abundance,colour=enrichment),alpha=.75)+
  scale_size_continuous(range=c(.12,7),guide="none")+
  scale_colour_manual(values=c(NotSig="grey80",Nodule="#d84c4c",Rhizosphere="#ff8d17",Root="#449f72"),breaks=c("Nodule","Rhizosphere","Root"),name="Enriched OTUs")+
  labs(x="Root",y="Nodule",z="Rhizosphere")+theme_bw(base_size=9)+theme(tern.axis.text=element_blank(),tern.axis.ticks=element_blank(),legend.position="right")
save_plot(p3,file.path(fig,"10c_microbe_niche_bubbles"),180,155)
## 4. 施肥条件：属水平颜色，绝非由图形推断丰度变化 ====
fert <- read.csv(file.path(dat,"fertilizer_otu.csv"))
genus_cols <- c(Klebsiella="#E31A1C",Pseudomonas="#228B22",Enterobacteriaceae="#1F78B4",Rosenbergiella="#FDB462",Oxalobacteraceae="#8B658B",Sphingobacterium="#4876FF",Lactococcus="#00BFFF",Erwinia="#EE82EE",Others="#8B8682",NotSig="#CDC9C9")
p4 <- ggtern(fert,aes(CK,NPK,NPKM))+geom_point(aes(size=abundance,colour=Genus),alpha=.8)+scale_size_continuous(range=c(.2,7),guide="none")+
  scale_colour_manual(values=genus_cols,guide=guide_legend(nrow=3,override.aes=list(size=2)))+theme_bw(base_size=9)+theme(legend.position="bottom",legend.text=element_text(size=6))
save_plot(p4,file.path(fig,"10d_fertilizer_genus_bubbles"),175,175)
## 5. 状态概率密度：投影后二维高斯核密度，三角外严格掩膜 ====
prob <- read.csv(file.path(dat,"state_probabilities.csv"))
pol <- data.frame(x=c(0,.5,1,0),y=c(0,sqrt(3)/2,0,0))
for(cond in c("OKN1elf","OKG148Le")) {
  dd <- prob[prob$condition==cond,]; xx <- dd$Noncanonical+.5*dd$Canonical; yy <- sqrt(3)/2*dd$Canonical
  kde <- MASS::kde2d(xx,yy,n=160,lims=c(0,1,0,sqrt(3)/2))
  zz <- expand.grid(x=kde$x,y=kde$y); zz$density <- as.vector(kde$z)
  zz <- zz[zz$x>=zz$y/sqrt(3) & zz$x<=1-zz$y/sqrt(3),]
  # 仅裁去密度极低尾部用于显示；不是论文5-95%概率质量等高线。
  zz$density[zz$density<max(zz$density)*.035] <- NA
  pal <- if(cond=="OKN1elf") c("white","#f5e9c8","#d8b365","#8c6d31") else c("white","#d5cfe5","#8c7bb5","#4b3c75")
  p5 <- ggplot(zz,aes(x,y))+geom_raster(aes(fill=density),na.rm=TRUE)+
    geom_path(data=pol,colour="grey40",linewidth=.4)+
    annotate("text",x=c(-.025,.5,1.04),y=c(-.035,.92,-.035),label=c("Stem-like","Canonical","Non-canonical"),size=2.3,colour=c("#329b68","grey20","#d75c68"))+
    scale_fill_gradientn(colours=pal,na.value="white",name="Density")+
    coord_fixed(xlim=c(-.12,1.12),ylim=c(-.10,.98))+labs(title=cond)+theme_void(base_size=9)+theme(plot.title=element_text(hjust=.5))
  save_plot(p5,file.path(fig,paste0("10e_state_density_",cond)),155,145)
}

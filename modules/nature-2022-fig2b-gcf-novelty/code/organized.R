# 01 独立 Fig.2b：官方 GCF 数据、图例和输出，不调用 a/c 或历史草稿 ----------
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
b <- "gallery/nature-2022-fig2b-gcf-novelty"
obs<-read.csv(file.path(b,"data/histogram-observations.csv"))
gcf<-read.csv(file.path(b,"data/gcf-distance-summary.csv"))
phyla<-c(Actinobacteriota="#F6BD82",Proteobacteria="#B3C7E5",Firmicutes="#C2B2D3",Cyanobacteria="#A8DC93",Bacteroidota="#F29D99",Thermoplasmatota="#BE9E96",Marinisomatota="#EEB9D1",Chloroflexota="#DCDA96",Verrucomicrobiota="#AAD9E3",Planctomycetota="#F9E8B9",`Other phyla`="#7F7F7F")
classes<-c(NRPS="#B3D568",T1PKS="#8DA0CB",`T2/3PKS`="#D68BC1",RiPPs="#EF8A62",Terpenes="#80BDA6",Other="#FED94B")
theme_set(theme_classic(base_size=8,base_family="Arial"))
# 02 与作者同口径：非零 GCF×category 各计一次；不把组成 weight 用作频数 ------
stopifnot(nrow(gcf)==6907,sum(gcf$RefSeq>.2)==3861,sum(gcf$MIBiG>.2)==6688)
obs$type<-factor(obs$type,levels=c(names(phyla),names(classes)))
obs$db<-factor(obs$db,levels=c("RefSeq","MIBiG"));obs$facet<-factor(obs$facet,levels=c("BGC class","Phyla"))
ann<-data.frame(facet=factor("BGC class",levels=levels(obs$facet)),db=factor("RefSeq",levels=levels(obs$db)))
lab_db<-setNames(paste0(c("RefSeq","MIBiG"),"\n",scales::comma(c(sum(gcf$RefSeq>.2),sum(gcf$MIBiG>.2)))," (",round(100*c(mean(gcf$RefSeq>.2),mean(gcf$MIBiG>.2))),"%)"),levels(obs$db))
p_hist<-ggplot(obs,aes(distance,fill=type))+geom_histogram(bins=30,colour="white",linewidth=.10)+geom_vline(xintercept=.2,linewidth=.3)+
 geom_segment(data=ann,aes(x=.22,xend=.96,y=1040,yend=1040),inherit.aes=FALSE,arrow=grid::arrow(ends="both",length=unit(1.1,"mm")),linewidth=.25)+
 geom_text(data=ann,aes(x=.59,y=1150,label="Novel (mean d > 0.2)"),inherit.aes=FALSE,size=2.6)+
 facet_grid(facet~db,labeller=labeller(db=lab_db))+
 scale_fill_manual(values=c(phyla,classes),drop=FALSE)+scale_x_continuous(breaks=c(0,.25,.5,.75,1),labels=c("0","0.25","0.50","0.75","1.00"),expand=expansion(mult=0))+
 scale_y_continuous(breaks=c(0,400,800,1200),labels=scales::label_comma(),expand=expansion(mult=0))+
 coord_cartesian(xlim=c(-.055,1.055),ylim=c(0,1300))+
 labs(x=expression("Mean minimum cosine distance of BGC members ("*bar(d)*")"),y="GCF-category count",tag="b")+
 theme(legend.position="none",strip.background=element_blank(),strip.text=element_text(size=8),strip.text.y=element_text(angle=-90),axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.25),axis.ticks.length=unit(.8,"mm"),axis.title=element_text(size=8),axis.text=element_text(size=7),panel.spacing.x=unit(6,"mm"),panel.spacing.y=unit(4,"mm"),plot.tag=element_text(face="bold",size=12),plot.tag.position=c(.005,.99),plot.margin=margin(2,3,2,2,"mm"))
# 03 本 panel 自带类别和门图例 -------------------------------------------
key<-bind_rows(data.frame(label=names(classes),fill=unname(classes),x=0,y=16-(0:5)),data.frame(label=names(phyla),fill=unname(phyla),x=0,y=8.7-.75*(0:10)))
p_leg<-ggplot()+geom_tile(data=key,aes(x=.2,y=y,fill=fill),width=.45,height=.45)+scale_fill_identity()+
 geom_text(data=key,aes(x=.6,y=y,label=label),size=2.45,hjust=0,family="Arial")+
 annotate("text",x=0,y=17,label="BGC class",hjust=0,size=2.8)+annotate("text",x=0,y=9.7,label="Phylum",hjust=0,size=2.8)+
 coord_cartesian(xlim=c(0,4.1),ylim=c(0,18),clip="off",expand=FALSE)+theme_void()+theme(plot.margin=margin(2,5,2,2,"mm"))
p<-p_hist+p_leg+plot_layout(widths=c(3.1,1))
p<-p+plot_annotation(caption="6,907 independent GCFs. Headline numbers are unique GCF counts; stacks count GCF-category assignments (author code).",theme=theme(plot.caption=element_text(size=6,hjust=0)))
# 04 输出与硬检查 -------------------------------------------------------
ggsave(file.path(b,"preview.png"),plot=p,device=ragg::agg_png,width=190,height=105,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"output/figures/fig2b-gcf-novelty.pdf"),plot=p,device=cairo_pdf,width=190,height=105,units="mm",bg="white")
ggsave(file.path(b,"validation/histograms.png"),plot=p_hist,device=ragg::agg_png,width=150,height=90,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"validation/legend.png"),plot=p_leg,device=ragg::agg_png,width=48,height=100,units="mm",dpi=300,bg="white")
build<-ggplot_build(p_hist);bins<-build$data[[1]]
stopifnot(sum(bins$count)==nrow(obs),!anyNA(obs$type),as.character(build$layout$layout$db[1])=="RefSeq",file.info(file.path(b,"preview.png"))$size>0)
write.csv(bins[c("PANEL","group","x","xmin","xmax","count","ymin","ymax","fill")],file.path(b,"data/histogram-bins.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(b,"validation/render-session.txt"))

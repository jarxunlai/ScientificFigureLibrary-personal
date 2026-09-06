# 01 独立 Fig.2c：真实 GTDB 拓扑 + BGC 堆叠径向柱 --------------------------
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
b <- "gallery/nature-2022-fig2c-phylogenomic-bgc"
seg <- read_csv(file.path(b,"data/tree-segments.csv"),show_col_types=FALSE)
wed <- read_csv(file.path(b,"data/tree-wedges.csv"),show_col_types=FALSE)
bins <- read_csv(file.path(b,"data/bgc-radial-bins.csv"),show_col_types=FALSE)
lab <- read_csv(file.path(b,"data/taxonomy-labels.csv"),show_col_types=FALSE)
unplaced <- read_csv(file.path(b,"data/unplaced-species.csv"),show_col_types=FALSE)
theme_set(theme_void(base_size=8,base_family="Arial"))
class_cols <- c(RiPPs="#EF8A62",NRPS="#B3D568",T1PKS="#8DA0CB",T2_3PKS="#D68BC1",Terpenes="#80BDA6",Other="#FED94B")
class_names <- c("RiPPs","NRPS","T1PKS","T2/3PKS","Terpenes","Other")

# 02 覆盖率色阶只用于可见终端枝/折叠支；纯参考支=灰色 ----------------------
seg <- seg %>% group_by(domain,edge) %>% mutate(terminal=abs(last(sqrt(x*x+y*y))-1)<1e-7) %>% ungroup()
max_genomes <- max(wed$n_genomes,seg$n_genomes[seg$terminal])
pal <- viridisLite::viridis(256,begin=.1,end=.85,direction=-1)
ramp <- scales::col_numeric(pal,domain=c(0,log(max_genomes)))
wed$color <- ramp(log(pmax(wed$n_genomes,1)));wed$color[wed$n_genomes==0] <- "#D4D4D4"
seg$color <- "#C8C8C8"
ix <- seg$terminal & seg$n_genomes>0
seg$color[ix] <- ramp(log(seg$n_genomes[ix]))
# 03 BGC 柱：同一选中基因组的类别/products × regions，然后真正堆叠 --------
base_radius <- 1.115
unit_radius <- .008
bar <- bins %>% filter(!is.na(species)) %>% mutate(bin_id=paste(domain,bin,sep="_")) %>%
 pivot_longer(all_of(names(class_cols)),names_to="class",values_to="count") %>%
 mutate(class=factor(class,levels=rev(names(class_cols))),height=if_else(products>0,count/products*regions,0)) %>%
 arrange(bin_id,class) %>% group_by(bin_id) %>%
 mutate(upper=cumsum(height),lower=upper-height,r0=base_radius+unit_radius*lower,r1=base_radius+unit_radius*upper) %>% ungroup()
check <- bar %>% group_by(bin_id) %>% summarise(height=sum(height),regions=first(regions),.groups="drop")
stopifnot(all(abs(check$height-check$regions)<1e-9),!anyNA(bar$height))
polys <- list()
for(i in which(bar$height>0)) {
 z<-bar[i,];aa<-seq(z$theta-z$half_width,z$theta+z$half_width,length.out=4)
 polys[[length(polys)+1]]<-data.frame(x=c(z$r0*cos(aa),z$r1*cos(rev(aa))),y=c(z$r0*sin(aa),z$r1*sin(rev(aa))),id=i,color=unname(class_cols[as.character(z$class)]))
}
polys<-bind_rows(polys)
# 同一数值刻度用于细菌和古菌，不分别缩放。
rings <- expand.grid(value=c(0,5,10,15,20,30),theta=c(seq(0,75,length.out=200),seq(90,360,length.out=700))*pi/180)
rings$sector<-ifelse(rings$theta<=75*pi/180,"Archaea","Bacteria")
rings$x<-(base_radius+unit_radius*rings$value)*cos(rings$theta)
rings$y<-(base_radius+unit_radius*rings$value)*sin(rings$theta)
axis_labels<-data.frame(value=c(0,5,10,15,20,30),theta=82*pi/180)
axis_labels$x<-(base_radius+unit_radius*axis_labels$value)*cos(axis_labels$theta)
axis_labels$y<-(base_radius+unit_radius*axis_labels$value)*sin(axis_labels$theta)
# taxonomy 名称从真实树节点获得，沿树外缘切向书写。
lab$x<-1.048*cos(lab$theta);lab$y<-1.048*sin(lab$theta)
lab$angle<-((lab$theta*180/pi+90+90)%%180)-90

# 04 原图六类 BGC-rich 标记，连接实际 bin，不按截图伪造坐标 ---------------
patterns<-c("f__Sandaracinaceae","g__Tistrella;","p__Planctomycetota;","p__Eremiobacterota;","g__Rhodococcus;","g__Synechococcus")
highlights<-list()
for(i in seq_along(patterns)) {
 z<-bins %>% filter(!is.na(taxonomy),grepl(patterns[i],taxonomy),regions>15) %>% arrange(desc(regions)) %>% slice_head(n=1)
 stopifnot(nrow(z)==1)
 z$mark<-letters[i];z$color<-if(i==4)"#D7301F" else "#555555"
 z$r<-base_radius+unit_radius*z$regions
 highlights[[i]]<-z
}
highlights<-bind_rows(highlights)
tri<-list()
for(i in seq_len(nrow(highlights))) {
 z<-highlights[i,]; rr<-c(z$r+.02,z$r+.06,z$r+.06);aa<-c(z$theta,z$theta-.012,z$theta+.012)
 tri[[i]]<-data.frame(x=rr*cos(aa),y=rr*sin(aa),id=i,color=z$color)
}
tri<-bind_rows(tri)
highlights$x<-(highlights$r+.095)*cos(highlights$theta);highlights$y<-(highlights$r+.095)*sin(highlights$theta)
p_circle<-ggplot()+
 geom_path(data=rings,aes(x,y,group=interaction(value,sector)),linewidth=.15,colour="#888888")+
 geom_polygon(data=wed,aes(x,y,group=interaction(domain,node),fill=color),colour=NA)+
 geom_path(data=seg,aes(x,y,group=interaction(domain,edge),colour=color),linewidth=.12)+
 geom_polygon(data=polys,aes(x,y,group=id,fill=color),colour=NA)+
 geom_polygon(data=tri,aes(x,y,group=id,fill=color),colour=NA)+
 geom_text(data=highlights,aes(x,y,label=mark,colour=color),size=2.7,fontface="bold")+
 geom_text(data=lab,aes(x,y,label=label,angle=angle),size=2.35)+
 geom_text(data=axis_labels,aes(x,y,label=value),size=2.05)+
 annotate("text",x=-.32,y=-.20,label="Bacteria",size=4,angle=-45)+
 annotate("text",x=.38,y=.30,label="Archaea",size=4,angle=-45)+
 scale_fill_identity()+scale_colour_identity()+coord_equal(xlim=c(-1.43,1.43),ylim=c(-1.43,1.43),expand=FALSE)+
 labs(tag="c")+theme(plot.tag=element_text(size=12,face="bold"),plot.tag.position=c(.005,.99),plot.margin=margin(1,1,1,1,"mm"))
# 05 本 panel 的完整图例、参考灰色定义、缺失位置的诚实说明 ------------------
key<-data.frame(x=0,y=16-(0:5),color=unname(class_cols),label=class_names)
colors_key<-data.frame(x=seq(0,5,length.out=150),y=6.5,color=ramp(seq(0,log(max_genomes),length.out=150)))
nt<-c(1,10,100,1000)
number_key<-data.frame(x=5*log(nt)/log(max_genomes),y=5.7,label=scales::comma(nt))
rich_key<-data.frame(x=0,y=2.1-.65*(0:5),label=paste0(letters[1:6],"  ",c("Sandaracinaceae","Tistrella","Planctomycetota","Eremiobacterota","Rhodococcus","Synechococcus")))
p_leg<-ggplot()+geom_tile(data=key,aes(x=.15,y=y,fill=color),width=.4,height=.35)+
 geom_text(data=key,aes(x=.6,y=y,label=label),hjust=0,size=2.5)+
 annotate("text",x=0,y=17,label="BGC class",hjust=0,size=2.8)+
 annotate("text",x=0,y=9.3,label="Outer ring: BGC regions\n0, 5, 10, 15, 20, 30",hjust=0,size=2.4)+
 annotate("text",x=0,y=7.45,label="Number of OMD genomes",hjust=0,size=2.6)+
 geom_tile(data=colors_key,aes(x,y,fill=color),width=.04,height=.35)+geom_text(data=number_key,aes(x,y,label=label),size=2.15)+
 annotate("rect",xmin=0,xmax=.4,ymin=4.25,ymax=4.6,fill="#D4D4D4")+
 annotate("text",x=.6,y=4.4,label="Reference-only clades",hjust=0,size=2.35)+
 annotate("text",x=0,y=3.1,label="BGC-rich lineages (>15)",hjust=0,size=2.6)+
 geom_text(data=rich_key,aes(x,y,label=label),hjust=0,size=2.35)+
 annotate("text",x=0,y=-3.3,label=paste0("Unplaced in author tree:\n",paste0(unplaced$species,": ",unplaced$regions," BGCs",collapse="\n")),hjust=0,size=2.2)+
 scale_fill_identity()+coord_cartesian(xlim=c(0,5.8),ylim=c(-6,18),clip="off",expand=FALSE)+theme(plot.margin=margin(2,5,2,3,"mm"))
p<-p_circle+p_leg+plot_layout(widths=c(2.8,1))
p<-p+plot_annotation(caption="Author GTDB trees + OMD data; 7,574 mapped species clusters (33,489 genomes).\nFour later MAG clusters (10 genomes) lack tree placement. Collapse/layout reimplemented; no GTDB-Tk rerun.",theme=theme(plot.caption=element_text(size=6,hjust=0),plot.margin=margin(2,2,2,2,"mm")))
# 06 输出与验证：preview 由本脚本直接生成，禁止用历史图覆盖 -----------------
ggsave(file.path(b,"preview.png"),plot=p,device=ragg::agg_png,width=190,height=145,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"output/figures/fig2c-phylogenomic-bgc.pdf"),plot=p,device=cairo_pdf,width=190,height=145,units="mm",bg="white")
ggsave(file.path(b,"validation/circle.png"),plot=p_circle,device=ragg::agg_png,width=130,height=130,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"validation/legend.png"),plot=p_leg,device=ragg::agg_png,width=50,height=125,units="mm",dpi=300,bg="white")
write_csv(bar,file.path(b,"data/stacked-bgc-heights.csv"));write_csv(highlights,file.path(b,"data/highlight-positions.csv"))
write_csv(check,file.path(b,"validation/stack-height-checks.csv"))
stopifnot(file.info(file.path(b,"preview.png"))$size>0,all(wed$color[wed$n_genomes==0]=="#D4D4D4"))
writeLines(capture.output(sessionInfo()),file.path(b,"validation/render-session.txt"))

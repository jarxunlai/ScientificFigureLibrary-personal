# 海洋微生物生物合成潜力：Nature 2022 Fig.2a
# 输入来自官方 Supplementary Table 2 + recover-author-tree.py 恢复的作者聚类树。
# 线性脚本；从项目根目录运行。只写本条目的 data/output/validation。

# 01 加载绘图输入与风格 --------------------------------------------------
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(ape); library(patchwork)})
b <- "gallery/nature-2022-fig2a-gcc-overview"
gcc <- read.csv(file.path(b, "data/gcc-summary.csv"), check.names = FALSE)
tax <- read.csv(file.path(b, "data/gcc-phylum-composition.csv"))
cls <- read.csv(file.path(b, "data/gcc-class-composition.csv"))


phylum_colors <- c(Actinobacteriota="#F6BD82", Proteobacteria="#B3C7E5", Firmicutes="#C2B2D3", Cyanobacteria="#A8DC93", Bacteroidota="#F29D99", Thermoplasmatota="#BE9E96", Marinisomatota="#EEB9D1", Chloroflexota="#DCDA96", Verrucomicrobiota="#AAD9E3", Planctomycetota="#F9E8B9", `Other phyla`="#7F7F7F")
class_colors <- c(NRPS="#B3D568", T1PKS="#8DA0CB", `T2/3PKS`="#D68BC1", RiPPs="#EF8A62", Terpenes="#80BDA6", Other="#FED94B")
genome_colors <- c(`REFs/SAGs only`="#DFF3DB", `REFs/SAGs and MAGs`="#A8DCB5", `MAGs only`="#43A2CA")
colors <- c(class_colors, phylum_colors, genome_colors)
# figure-style 原则的 R 实现：三级字号、白底、显式设备、图例无边框。
theme_set(theme_classic(base_size=8, base_family="Arial"))

# 02 将作者聚类树修剪为 151 GCC，并在 GCC 定义阈值处终止枝条 ------------
tr <- read.tree(file.path(b, "data/author-gcc-151.newick"))
stopifnot(setequal(tr$tip.label, gcc$gcc), Ntip(tr)==151)
n <- Ntip(tr)
# 原模型内 GCC 子树的高度不同。只显示 GCC 以上的层次关系：所有 tip 在 d=0.8 处截断。
# 内部节点的 merge distance 和拓扑不变；这不是作者最终 rerooted Newick 的逐字恢复。
cluster_parameters <- read.csv(file.path(b, "data/author-clustering-parameters.csv"), stringsAsFactors=FALSE)
root_height <- as.numeric(cluster_parameters$value[cluster_parameters$parameter=="root_merge_distance"])
# 保留 151 tips 的根可能包含被剪去的 stem；用与原 197 树的共享 tip 深度差校正。
tr_full <- read.tree(file.path(b, "data/author-gcc-197.newick"))
dep_full <- node.depth.edgelength(tr_full)
dep <- node.depth.edgelength(tr)
shared <- tr$tip.label[1]
stem <- dep_full[match(shared, tr_full$tip.label)] - dep[1]
# 完整原始聚类树的 root merge distance 来自恢复脚本；根据作者 d=0.8 cut 指定叶端。
cut_depth <- root_height - stem - as.numeric(cluster_parameters$value[cluster_parameters$parameter=="gcc_threshold"])
for (i in which(tr$edge[,2] <= n)) {
  tr$edge.length[i] <- cut_depth - dep[tr$edge[i,1]]
}
stopifnot(all(tr$edge.length >= -1e-8))
tr$edge.length <- pmax(tr$edge.length,0)
tr <- ladderize(tr, right=TRUE)
tr <- reorder.phylo(tr,"cladewise")
write.tree(tr,file.path(b,"data/author-gcc-151-cut-0.8.newick"))
# 仅调换姊妹枝绘图顺序，不重新聚类、不人工改变数据归属。
leaf_order <- tr$edge[tr$edge[,2] <= n,2]
angle <- rep(NA_real_,n+tr$Nnode)
angle[leaf_order] <- pi/2 - (2.5+(seq_len(n)-0.5)*355/n)*pi/180
for (p in rev(unique(tr$edge[,1]))) angle[p] <- mean(angle[tr$edge[tr$edge[,1]==p,2]])
dep <- node.depth.edgelength(tr)
radius <- 0.12+0.88*dep/max(dep)
stopifnot(!anyNA(angle), max(radius)<=1.00001)
leaves <- data.frame(gcc=tr$tip.label,theta=angle[seq_len(n)], tip=n:1)
gcc <- left_join(gcc,leaves,by="gcc")
stopifnot(!anyNA(gcc$theta))
edges <- list()
for(i in seq_len(nrow(tr$edge))) {
  p <- tr$edge[i,1]; ch <- tr$edge[i,2]
  aa <- seq(angle[p],angle[ch],length.out=18)
  edges[[i]] <- data.frame(x=c(radius[p]*cos(aa),radius[ch]*cos(angle[ch])), y=c(radius[p]*sin(aa),radius[ch]*sin(angle[ch])), edge=i)
}
edges <- bind_rows(edges)

# 03 准备圆环的真实数值与 ID 一一对应关系 --------------------------------
red_ramp <- scales::col_numeric(c("#FEE5D9","#CB181D"),domain=c(0,.65))
blue_ramp <- scales::col_numeric(c("#ECE7F2","#0570B0"),domain=c(0,1))
grey_ramp <- scales::col_numeric(c("#EEEEEE","#222222"),domain=c(0,log1p(max(gcc$n_bgc))))
track <- bind_rows(
  transmute(gcc,gcc,theta,r0=1.045,r1=1.082,fill=unname(genome_colors[genome_type])),
  transmute(gcc,gcc,theta,r0=1.955,r1=2.005,fill=grey_ramp(log1p(n_bgc))),
  transmute(gcc,gcc,theta,r0=2.018,r1=2.103,fill=blue_ramp(prevalence)),
  transmute(gcc,gcc,theta,r0=2.115,r1=2.240,fill=red_ramp(refseq)))
tax <- tax %>% mutate(phylum_plot=factor(phylum_plot,levels=rev(names(phylum_colors)))) %>% arrange(gcc,phylum_plot) %>% group_by(gcc) %>% mutate(upper=cumsum(prop),lower=upper-prop) %>% ungroup() %>% left_join(leaves,by="gcc")
stopifnot(all(abs(tapply(tax$prop,tax$gcc,sum)-1)<1e-8))
track <- bind_rows(track,transmute(tax,gcc,theta,r0=1.115+.53*lower,r1=1.115+.53*upper,fill=unname(phylum_colors[as.character(phylum_plot)])))
# taxonomy 扇形较窄；其余 heatmap track 近乎相接。
track$half <- ifelse(track$r0>=1.11 & track$r1<=1.65,.36,.49)*355/n*pi/180
polygons <- list()
for(i in seq_len(nrow(track))) {
  z <- track[i,]; aa <- seq(z$theta-z$half,z$theta+z$half,length.out=4)
  polygons[[i]] <- data.frame(x=c(z$r0*cos(aa),z$r1*cos(rev(aa))),y=c(z$r0*sin(aa),z$r1*sin(rev(aa))),fill=z$fill,id=i)
}
polygons <- bind_rows(polygons)
cls <- left_join(cls,leaves,by="gcc")
cls$radius <- 1.695+.044*(match(cls$class,names(class_colors))-1)
cls$x <- cls$radius*cos(cls$theta); cls$y <- cls$radius*sin(cls$theta)
stopifnot(!anyNA(cls$x), !anyNA(polygons$fill), all(gcc$prevalence>=0 & gcc$prevalence<=1))
# 参考线只表示类别 track，点面积表示 GCC 内该类别频率。
guide_circles <- expand.grid(r=1.695+.044*(0:5),a=seq(92.5,447.5,length.out=600)*pi/180)
guide_circles$x <- guide_circles$r*cos(guide_circles$a); guide_circles$y <- guide_circles$r*sin(guide_circles$a)
# 最小 MIBiG 距离：与图例一致的黑/灰箭头；显式对应真实 GCC。
hits <- gcc[gcc$mibig_min<.1,]
triangles <- list()
for(i in seq_len(nrow(hits))) {
  a <- hits$theta[i]; rr <- c(2.272,2.329,2.329); aa <- c(a,a-.015,a+.015)
  triangles[[i]] <- data.frame(x=rr*cos(aa),y=rr*sin(aa),fill=ifelse(hits$mibig_min[i]<1e-10,"black","#7F7F7F"),id=i)
}
triangles <- bind_rows(triangles)
# Nature 原图的 A-D 标签来自手工标注的 GCC/天然产物映射；本输入表没有该字段。
# 只把能由作者脚本注释直接核验的映射写入；不把大型 GCC 的距离最小值误当成 arylpolyene。
known <- data.frame(gcc=c("gcc_92","gcc_123","gcc_14","gcc_22"),label=c("A","B","C","D")) %>% left_join(leaves,by="gcc")
known$x <- 2.40*cos(known$theta); known$y <- 2.40*sin(known$theta)
# A: GCC00092 arylpolyene；B: GCC00123 ectoine；C: GCC00014 carotenoid；D: GCC00022 siderophore。
track_labels <- data.frame(y=c(2.18,2.06,1.82,1.39,1.064),label=as.character(1:5))
p_a <- ggplot()+
  geom_path(data=guide_circles,aes(x,y,group=r),colour="grey80",linewidth=.11)+
  geom_polygon(data=polygons,aes(x,y,group=id,fill=fill),colour="white",linewidth=.06)+
  geom_path(data=edges,aes(x,y,group=edge),linewidth=.12,colour="black")+
  geom_point(data=cls,aes(x,y,size=prop,fill=unname(class_colors[class])),shape=21,stroke=0)+
  scale_size_area(max_size=1.20,limits=c(0,1),guide="none")+
  geom_polygon(data=triangles,aes(x,y,group=id,fill=fill),colour=NA)+
  geom_text(data=known,aes(x,y,label=label),size=2.65,family="Arial")+
  geom_text(data=track_labels,aes(x=0,y=y,label=label),size=2.6,family="Arial")+
  scale_fill_identity()+coord_equal(xlim=c(-2.47,2.47),ylim=c(-2.46,2.49),expand=FALSE)+theme_void()+
  labs(tag="a")+theme(plot.tag=element_text(face="bold",size=12),plot.tag.position=c(.015,.99),plot.margin=margin(1,1,1,1,"mm"))

# 04 信息图例：色条、分类色块、圆点大小、来源类型及 MIBiG 标记 ---------------
legend_text <- data.frame(x=numeric(),y=numeric(),label=character(),size=numeric())
legend_text <- bind_rows(legend_text,
 data.frame(x=2,y=c(124,109,94,77,25),label=c("1  GCC novelty\n    Mean RefSeq distance","2  Prevalence (%)","BGC count","3  BGC class","Distance to known BGCs"),size=8),
 data.frame(x=57,y=c(124,117,87,35),label=c("4  BGC distribution (phylum)","Frequent in MIBiG (>5%)","Frequent in OMD (>1%)","5  Genome types"),size=c(8,7,7,8)))
swatches <- bind_rows(
 data.frame(x=57,y=110-6*(0:3),label=names(phylum_colors)[1:4],fill=unname(phylum_colors[1:4])),
 data.frame(x=57,y=80-6*(0:6),label=c(names(phylum_colors)[5:10],"Other"),fill=unname(phylum_colors[5:11])),
 data.frame(x=57,y=27-7*(0:2),label=names(genome_colors),fill=unname(genome_colors)))
legend_text <- bind_rows(legend_text,transmute(swatches,x=x+5,y,label,size=7))
class_key <- data.frame(x=3,y=70-6*(0:5),label=names(class_colors),fill=unname(class_colors))
legend_text <- bind_rows(legend_text,transmute(class_key,x=x+6,y,label,size=7))
# 分频泡泡图例额外补充原图未量化的 size 标尺。
bubble_key <- data.frame(x=c(7,20,33),y=33,p=c(.1,.5,1),label=c("10%","50%","100%"))
legend_text <- bind_rows(legend_text,transmute(bubble_key,x=x-3,y=y-4,label,size=6),
 data.frame(x=c(7,7,2,2),y=c(21,16,9,3),label=c("min(MIBiG d) < 1e-10","min(MIBiG d) < 0.1","A  Arylpolyene   B  Ectoine","C  Carotenoid    D  Siderophore"),size=6.5))
gradients <- bind_rows(data.frame(x=seq(3,44,length.out=100),y=117,fill=red_ramp(seq(0,.65,length.out=100))),data.frame(x=seq(3,44,length.out=100),y=102,fill=blue_ramp(seq(0,1,length.out=100))),data.frame(x=seq(3,44,length.out=100),y=87,fill=grey_ramp(seq(0,log1p(max(gcc$n_bgc)),length.out=100))))
legend_text <- bind_rows(legend_text,
 data.frame(x=3+41*c(0,.2,.4,.6)/.65,y=112.5,label=c("0","0.2","0.4","0.6"),size=6),
 data.frame(x=3+41*c(0,.25,.5,.75,1),y=97.5,label=c("0","25","50","75","100"),size=6),
 data.frame(x=3+41*log1p(c(1,50,400,3000))/log1p(max(gcc$n_bgc)),y=82.5,label=c("1","50","400","3,000"),size=6))
p_leg <- ggplot()+geom_tile(data=gradients,aes(x,y,fill=fill),width=.5,height=2.7)+
 geom_tile(data=swatches,aes(x=x+1,y=y,fill=fill),width=2.8,height=3.0)+
 geom_point(data=class_key,aes(x,y,fill=fill),shape=21,size=2.5,stroke=0)+
 geom_point(data=bubble_key,aes(x,y,size=p),colour="#555555")+scale_size_area(max_size=1.20,limits=c(0,1),guide="none")+
 annotate("point",x=c(3,3),y=c(21,16),shape=17,size=2.1,colour=c("black","#7F7F7F"))+
 geom_text(data=legend_text,aes(x,y,label=label,size=I(size/2.845)),hjust=0,family="Arial")+scale_fill_identity()+
 coord_cartesian(xlim=c(0,112),ylim=c(-2,127),clip="off",expand=FALSE)+theme_void()+theme(plot.margin=margin(2,5,2,2,"mm"))

# 05 只输出 panel a；不再读写 panel b 或旧组合草稿 -----------------------
p <- p_a+p_leg+plot_layout(widths=c(1.13,1))
p <- p+plot_annotation(caption="151 GCCs; 17,689 representative BGCs. Author clustering recovered; final reroot/leaf order unavailable.",theme=theme(plot.caption=element_text(size=6,hjust=0)))
ggsave(file.path(b,"preview.png"),plot=p,device=ragg::agg_png,width=190,height=117,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"output/figures/fig2a-gcc-overview.pdf"),plot=p,device=cairo_pdf,width=190,height=117,units="mm",bg="white")
ggsave(file.path(b,"validation/circle.png"),plot=p_a,device=ragg::agg_png,width=110,height=110,units="mm",dpi=300,bg="white")
ggsave(file.path(b,"validation/legend.png"),plot=p_leg,device=ragg::agg_png,width=90,height=112,units="mm",dpi=300,bg="white")
write.csv(gcc[c("gcc","n_bgc","prevalence","refseq","mibig_min","genome_type","theta")],file.path(b,"data/position-map.csv"),row.names=FALSE)
checks<-data.frame(metric=c("GCC","representative_BGC","MAG_only_GCC","novel_GCC_RefSeq_gt_0.4"),observed=c(nrow(gcc),sum(gcc$n_bgc),sum(gcc$genome_type=="MAGs only"),sum(gcc$refseq>.4)),expected=c(151,17689,53,44))
stopifnot(all(checks$observed==checks$expected),file.info(file.path(b,"preview.png"))$size>0)
write.csv(checks,file.path(b,"validation/data-checks.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(b,"validation/render-session.txt"))
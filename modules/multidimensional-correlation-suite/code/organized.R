# 一维二维多维相关性可视化 / One-, two- and multidimensional correlation visualizations
# 科学问题：单变量关联、跨矩阵关联及群落-环境距离关联如何分别展示？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "multidimensional-correlation-suite")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：去除来源中没有的强制3×3切块，恢复连续双聚类热图、注释和Mantel P区间语义；输入/派生值均为合成示例，不冒充来源结果。
# 第 7 篇包含三个独立完整变体，全部执行；不把三种数据当作同一个分析。
# a: log10 观测的 Pearson 检验；b: 明确的示例跨矩阵相关值；c: 模拟样本的 partial Mantel。
# 作者代码中对相关矩阵再做 cor.mtest 的问题不沿用；检验基于原始观测的 r 与实际 n。
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(vegan)


source(file.path(module_root,"code/export.R"))

## 1. 一维变体：分组变量相关矩阵，P>=0.01 留空但保留格线 ====
obs<-read.csv(file.path(module_root,"data/environment_observations.csv"),check.names=FALSE)
grp<-read.csv(file.path(module_root,"data/environment_groups.csv"))
x<-log10(as.matrix(obs[,-1])+1e-6)
r<-cor(x); n<-ncol(r); df<-nrow(x)-2
p<-2*pt(-abs(r*sqrt(df/pmax(1-r^2,1e-15))),df)
diag(p)<-0
col_fun<-colorRamp2(c(-1,0,1),c("#b2182b","white","#1a1a1a"))
gcol<-c(Microbial="#66009b",Site="#00bccc","Land use"="#ff9933",Watershed="#00089d")
write.csv(r,file.path(draft,"result/07a_pearson_r.csv"));write.csv(p,file.path(draft,"result/07a_pearson_p.csv"))
save_draw({
  par(mar=c(.1,.1,.1,.1),family="Arial",xaxs="i",yaxs="i")
  plot.new();plot.window(xlim=c(-16,n+22),ylim=c(-7,n+19),asp=1)
  for(i in 1:n) for(j in i:n) {
    yy<-n-i+1
    rect(j-.5,yy-.5,j+.5,yy+.5,col="white",border="#c6c6c6",lwd=.3)
    if(p[i,j]<.01) {
      half<-.46*sqrt(abs(r[i,j]))
      rect(j-half,yy-half,j+half,yy+half,col=col_fun(r[i,j]),border=NA)
    }
  }
  for(i in 1:n) {
    text(i-.65,n-i+1,colnames(r)[i],adj=1,cex=.56,col=gcol[grp$Group[i]])
    text(i,n+.8,colnames(r)[i],srt=90,adj=0,cex=.57,col=gcol[grp$Group[i]])
    rect(n+1,n-i+.5,n+2.1,n-i+1.5,col=gcol[grp$Group[i]],border=NA)
    text(n+2.7,n-i+1,colnames(r)[i],adj=0,cex=.56)
  }
  for(g in unique(grp$Group)) {
    ids<-which(grp$Group==g)
    text(n+1.55,n-mean(ids)+1,g,srt=90,col=ifelse(g%in%c("Microbial","Watershed"),"white","black"),cex=.46)
  }
  z<-seq(-1,1,length.out=101)
  rect(8+head(z+1,-1)*10,-3.5,8+tail(z+1,-1)*10,-2.4,col=col_fun(head(z,-1)),border=NA)
  text(8+(c(-1,-.5,0,.5,1)+1)*10,-4.5,c(-1,-.5,0,.5,1),cex=.5)
  text(18,-6.1,"Pearson r; blank: P >= 0.01",cex=.52)
},file.path(draft,"figure/07a_grouped_correlation"),180,178)

## 2. 二维变体：跨矩阵相关值 + 双聚类树 + 细胞谱系注释 ====
cross<-as.matrix(read.csv(file.path(module_root,"data/cell_factor_correlations.csv"),row.names=1,check.names=FALSE))
cross<-cross*.55 # 合成跨矩阵减弱幅度，保存同一绘图输入；不是来源论文的相关值。
cg<-read.csv(file.path(module_root,"data/cell_groups.csv"))
cg_col<-c(Endothelium="#b36692",Epithelium="#fbb03c",Lymphoid="#0070bc",Mesenchyme="#faee90",Myeloid="#a0e0b8")
# scico 未安装，直接保留已观察到的青绿-白-粉紫色阶，不安装额外依赖。
cross_col<-colorRamp2(c(-.8,-.4,-.05,0,.05,.4,.8),c("#123d35","#72b5a1","#e4f2ed","white","#f2e7ef","#bc8ab8","#753a75"))
ha<-rowAnnotation(cell_group=cg$cell_group,col=list(cell_group=cg_col),show_annotation_name=FALSE,
                  simple_anno_size=unit(3,"mm"),gp=gpar(col="#999999",lwd=.3),
                  annotation_legend_param=list(labels_gp=gpar(fontsize=8),title_gp=gpar(fontsize=8)))
ht<-Heatmap(cross,name="r",col=cross_col,left_annotation=ha,cluster_rows=TRUE,cluster_columns=TRUE,

             row_title=NULL,column_title=NULL,row_names_gp=gpar(fontsize=7.5,fontfamily="Helvetica"),
             column_names_gp=gpar(fontsize=7.5,fontfamily="Helvetica"),row_names_side="right",
             rect_gp=gpar(col="#aaaaaa",lwd=.35),row_dend_width=unit(6,"mm"),column_dend_height=unit(6,"mm"),
             heatmap_legend_param=list(at=c(-.8,-.4,0,.4,.8),title_gp=gpar(fontsize=8),labels_gp=gpar(fontsize=8)),use_raster=FALSE)
set.seed(20261026)
save_draw({set.seed(20261026);draw(ht,padding=unit(c(3,3,3,3),"mm"),merge_legends=TRUE)},file.path(draft,"figure/07b_cross_matrix_heatmap"),180,143)
write.csv(cross,file.path(draft,"result/07b_plot_input.csv"))

## 3. 多维变体：由独立模拟样本实际计算 partial Mantel（999 次置换） ====
set.seed(20261026)
env<-read.csv(file.path(module_root,"data/mantel_environment.csv"),check.names=FALSE)
geo<-read.csv(file.path(module_root,"data/geography.csv"))
# 示例经纬度以按列标准化的欧氏距离表示，仅用于展示控制距离矩阵流程。
geo_dist<-dist(scale(geo[,c("Longitude","Latitude")]))
mantel<-data.frame()
community_inputs<-list(
 host=read.csv(file.path(module_root,"data/host.csv")),
 mOTUs=read.csv(file.path(module_root,"data/mOTUs.csv")),
 virus=read.csv(file.path(module_root,"data/virus.csv")),
 vOTUs=read.csv(file.path(module_root,"data/vOTUs.csv")))
for(g in c("host","mOTUs","virus","vOTUs")) {
  community<-community_inputs[[g]]
  community_dist<-vegdist(community[,-1],method="bray")
  for(v in names(env)[-1]) {
    test<-mantel.partial(community_dist,dist(scale(env[,v])),geo_dist,method="pearson",permutations=999)
    mantel<-rbind(mantel,data.frame(Community=g,Variable=v,r=unname(test$statistic),P=test$signif))
  }
}
write.csv(mantel,file.path(draft,"result/07c_partial_mantel.csv"),row.names=FALSE)
er<-cor(env[,-1]); en<-ncol(er); edf<-nrow(env)-2
ep<-2*pt(-abs(er*sqrt(edf/pmax(1-er^2,1e-15))),edf)
write.csv(er,file.path(draft,"result/07c_environment_r.csv"));write.csv(ep,file.path(draft,"result/07c_environment_p.csv"))
cells<-expand.grid(i=1:en,j=1:en); cells<-cells[cells$i<cells$j,]
cells$r<-er[cbind(cells$i,cells$j)];cells$P<-ep[cbind(cells$i,cells$j)]
cells$y<-en-cells$i+1
cells$half<-.43*sqrt(abs(cells$r))
cells$star<-ifelse(cells$P<.001,"***",ifelse(cells$P<.01,"**",ifelse(cells$P<.05,"*","")))
nodes<-data.frame(Community=c("host","mOTUs","virus","vOTUs"),x=c(-2.7,-.7,1.2,3.2),y=c(7.9,5.9,3.8,1.6))
mantel$x<-nodes$x[match(mantel$Community,nodes$Community)]
mantel$y<-nodes$y[match(mantel$Community,nodes$Community)]
mantel$xend<-match(mantel$Variable,names(env)[-1])
mantel$yend<-en-mantel$xend+1
mantel$pcat<-cut(mantel$P,c(-Inf,.001,.05,Inf),right=TRUE,labels=c("<= 0.001","0.001 - 0.05","> 0.05"))
mantel$rcat<-cut(mantel$r,c(-Inf,.1,.2,Inf),right=FALSE,labels=c("< 0.1","0.1 - 0.2",">= 0.2"))
endpoints<-data.frame(x=1:en,y=en:1,label=names(env)[-1])
p3<-ggplot()+
  geom_tile(data=cells,aes(j,y),fill="white",colour="#d0d0d0",linewidth=.25)+
  geom_rect(data=cells,aes(xmin=j-half,xmax=j+half,ymin=y-half,ymax=y+half,fill=r),colour=NA)+
  geom_text(data=cells,aes(j,y,label=star),size=2.4)+
  geom_curve(data=mantel,aes(x=x,y=y,xend=xend,yend=yend,colour=pcat,linewidth=rcat),curvature=.12,alpha=.8)+
  geom_point(data=endpoints,aes(x,y),shape=21,fill="#9ba9f1",colour="blue",size=1.5,stroke=.35)+
  geom_point(data=nodes,aes(x,y),shape=21,fill="#9ba9f1",colour="blue",size=1.6,stroke=.35)+
  geom_text(data=nodes,aes(x-.2,y,label=Community),hjust=1,size=2.8)+
  geom_text(data=endpoints,aes(x=12.1,y=y,label=label),hjust=0,size=2.8)+
  scale_x_continuous(breaks=1:en,labels=names(env)[-1],position="top",expand=c(0,0))+
  scale_y_continuous(breaks=NULL,expand=c(0,0))+
  scale_fill_gradient2(low="#FF8040",mid="white",high="#5BC2CD",limits=c(-1,1),name="Pearson r")+
  scale_colour_manual(values=c("#87CEEB","#9ACD32","#bfbfbf"),drop=FALSE,name="Mantel P")+
  scale_linewidth_manual(values=c(.18,.4,.7),drop=FALSE,name="Mantel r")+
  coord_fixed(xlim=c(-4.2,15.2),ylim=c(.25,11.75),clip="off")+
  labs(x=NULL,y=NULL)+theme_void(base_family="Arial",base_size=8)+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=.5,colour="black",size=8),
        legend.position="right",legend.text=element_text(size=7),legend.title=element_text(size=8),
        plot.margin=margin(5,3,2,2,"mm"))+
  guides(colour=guide_legend(order=1,override.aes=list(linewidth=.5)),linewidth=guide_legend(order=2),fill=guide_colourbar(order=3,barheight=unit(23,"mm")))
save_gg(p3,file.path(draft,"figure/07c_mantel_network"),180,125)

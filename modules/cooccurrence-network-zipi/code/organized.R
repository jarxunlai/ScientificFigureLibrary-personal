# 微生物共现网络与Zi-Pi角色 / Co-occurrence networks and Zi-Pi roles
# 科学问题：OTU相关网络中哪些节点跨模块连接较多，哪些节点模块内连接较强？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "cooccurrence-network-zipi")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 微生物网络与 Zi-Pi：恢复 A/B 网络 + C/D 角色图的原始四面板布局。
# 仅合成示例，保留现有计数表；为视觉示例生成更丰富的相关模块与桥接 OTU。
# 仍然从计数计算相关、筛边和角色；不手填 Zi/Pi，不把示例拓扑作为研究结论。
Sys.setlocale("LC_CTYPE","English_United States.utf8")



source(file.path(module_root,"code/export.R"))
options(warn=1)
options(device=function(...) grDevices::cairo_pdf(file=tempfile(fileext=".pdf"),family="Arial",...))
library(igraph)
fig <- file.path(entry,"figure"); res <- file.path(entry,"result")
phylum_cols <- c(Ascomycota="#5e4fa2",Basidiomycota="#ffff99",Unclassified="#b3a2d0",
  Mortierellomycota="#f2786d",Rozellomycota="#80b1d3",Chytridiomycota="#fdb462",
  Kickxellomycota="#b3de69",Zoopagomycota="#fccde5",Olpidiomycota="#d9d9d9")
networks <- list(); role_plots <- list(); metrics <- list()

## 1. 现有合成计数提取秩结构，生成模块/桥接混合的示例计数 ====
otu_inputs <- list(
  G=read.csv(file.path(module_root,"data/otu_counts_G.csv"),row.names=1),
  T=read.csv(file.path(module_root,"data/otu_counts_T.csv"),row.names=1))
for(group in c("G","T")) {
  set.seed(if(group=="G") 260113 else 260114)
  raw <- otu_inputs[[group]]
  tax <- read.csv(file.path(module_root,"data/taxonomy.csv"))
  n <- if(group=="G") 560 else 295
  # 节点扩展只服务示例网络密度；每一个新增节点、计数都明确另存 result。
  ns <- ncol(raw); module_n <- 6
  factors <- t(prcomp(scale(t(log1p(as.matrix(raw)))))$x[,2:7])
  factors <- t(scale(t(factors)))
  common <- as.numeric(scale(prcomp(scale(t(log1p(as.matrix(raw)))))$x[,1]))
  modules <- rep(1:module_n,length.out=n)
  counts <- matrix(0,n,ns)
  for(i in 1:n) {
    m <- modules[i]
    # 重叠载荷使网络不是六个断开的星状团；相关阈值仍保持作者的 |r|>.6。
    weights <- rep(0,module_n); weights[m] <- .85; weights[(m%%module_n)+1] <- .2
    if(i%%5==0) weights[] <- .16
    direction <- if(i%%23==0) -1 else 1
    latent <- direction*.9*common+as.numeric(weights %*% factors)+rnorm(ns,sd=runif(1,.35,.75))
    counts[i,] <- pmax(1,round(exp(8+.55*latent)))
  }
  rownames(counts) <- sprintf("SYN_%s_OTU%04d",group,seq_len(n))
  colnames(counts) <- colnames(raw)
  taxonomy <- data.frame(OTU_ID=rownames(counts),
    Phylum=sample(names(phylum_cols),n,replace=TRUE,prob=c(.43,.14,.13,.08,.06,.06,.04,.03,.03)))
  # 合成固定总量背景避免比例闭合把共同因子完全抵消；背景不作为真菌 OTU 节点。
  full_counts <- rbind(counts,SYN_BACKGROUND=max(colSums(counts))*6-colSums(counts))
  write.csv(full_counts,file.path(res,paste0("13_",group,"_visual_example_counts.csv")))
  write.csv(taxonomy,file.path(res,paste0("13_",group,"_visual_example_taxonomy.csv")),row.names=FALSE)

  ## 2. 真实计算此合成输入的边、模块与 Zi/Pi，不从来源图抄结果 ====
  prop <- sweep(counts,2,colSums(full_counts),"/")
  cr <- Hmisc::rcorr(t(sqrt(prop)),type="spearman")
  adj <- cr$r; adj[abs(cr$r)<.6 | cr$P>.05] <- 0; diag(adj) <- 0
  g <- graph_from_adjacency_matrix(adj,mode="undirected",weighted=TRUE,diag=FALSE)
  g <- delete_vertices(g,V(g)[degree(g)==0])
  E(g)$rho <- E(g)$weight; E(g)$weight <- abs(E(g)$weight)
  comm <- cluster_fast_greedy(g,weights=E(g)$weight); mem <- membership(comm)
  aa <- as_adjacency_matrix(g,attr=NULL,sparse=FALSE)
  kin <- sapply(sort(unique(mem)),function(m) rowSums(aa[,mem==m,drop=FALSE]))
  deg <- rowSums(aa); pi_value <- 1-rowSums((kin/deg)^2)
  own <- kin[cbind(seq_len(nrow(kin)),match(mem,sort(unique(mem))))]
  z <- numeric(vcount(g))
  for(m in sort(unique(mem))) {
    ids <- which(mem==m); sd_own <- sd(own[ids])
    z[ids] <- if(length(ids)>1 && is.finite(sd_own) && sd_own>0) (own[ids]-mean(own[ids]))/sd_own else NA_real_
  }
  stats <- data.frame(OTU_ID=V(g)$name,Degree=deg,Module=as.numeric(mem),Pi=pi_value,Zi=z)
  stats$Role <- ifelse(is.na(z),"Undefined Zi",ifelse(z>2.5,ifelse(pi_value>.62,"Network hubs","Module hubs"),ifelse(pi_value>.62,"Connectors","Peripherals")))
  write.csv(stats,file.path(res,paste0("13_",group,"_node_roles.csv")),row.names=FALSE)
  write_graph(g,file.path(res,paste0("13_",group,"_network.graphml")),format="graphml")
  metrics[[group]] <- data.frame(group=group,nodes=vcount(g),edges=ecount(g),
    positive=sum(E(g)$rho>0),negative=sum(E(g)$rho<0),modularity=modularity(comm),
    connectors=sum(stats$Role=="Connectors"),module_hubs=sum(stats$Role=="Module hubs"),undefined_zi=sum(is.na(z)))

  ## 3. 恢复紧凑球状网络、弧边和按门着色的可见节点 ====
  coords <- layout_with_fr(g,weights=E(g)$weight,niter=1000,grid="nogrid")
  coords <- scale(coords,center=TRUE,scale=FALSE)
  radii <- sqrt(rowSums(coords^2)); angles <- atan2(coords[,2],coords[,1])
  # FR 保留角度邻接顺序，再用分位数展开为紧凑圆盘，避免被几根长边拉成星形。
  # 这是可视化布局选择，不是对节点拓扑/统计结果的改写。
  theta <- 2*pi*(rank(angles,ties.method="first")-.5)/length(angles)-pi
  new_r <- sqrt((rank(radii,ties.method="first")-.5)/length(radii))*.98
  coords <- cbind(cos(theta)*new_r,sin(theta)*new_r)
  el <- as_edgelist(g,names=FALSE)
  edges <- data.frame(x=coords[el[,1],1],y=coords[el[,1],2],
    xend=coords[el[,2],1],yend=coords[el[,2],2],sign=ifelse(E(g)$rho>0,"Positive","Negative"))
  nodes <- data.frame(x=coords[,1],y=coords[,2],Phylum=taxonomy$Phylum[match(V(g)$name,taxonomy$OTU_ID)],Degree=deg)
  write.csv(nodes,file.path(res,paste0("13_",group,"_network_layout.csv")),row.names=FALSE)
  pn <- ggplot()+geom_curve(data=edges,aes(x,y,xend=xend,yend=yend,colour=sign),
      curvature=.12,linewidth=.12,alpha=if(group=="G") .055 else .15)+
    scale_colour_manual(values=c(Positive="#f18474",Negative="#70b79d"),guide="none")+
    ggnewscale::new_scale_colour()+
    geom_point(data=nodes,aes(x,y,size=Degree,fill=Phylum),shape=21,colour="#777777",stroke=.12)+
    scale_fill_manual(values=phylum_cols,breaks=names(phylum_cols),drop=FALSE,
      labels=c("Ascomycota","Basidiomycota","unclassified_k__Fungi","Mortierellomycota","Rozellomycota",
               "Chytridiomycota","Kickxellomycota","Zoopagomycota","Olpidiomycota"))+
    scale_size_continuous(range=c(.45,2.25),guide="none")+
    coord_equal(xlim=c(-1.05,1.05),ylim=c(-1.05,1.05),expand=FALSE)+
    theme_void(base_size=10,base_family="serif")+
    theme(legend.title=element_blank(),legend.key.height=unit(4.3,"mm"),legend.key.width=unit(4,"mm"),
      plot.margin=margin(1,1,1,1))+
    guides(fill=guide_legend(override.aes=list(size=3)))
  networks[[group]] <- pn
  save_plot(pn,file.path(fig,paste0("13_",group,"_cooccurrence_network")),155,135)

  ## 4. 角色图阈值与象限文字恢复；形状与作者图一致 ====
  stats$Role <- factor(stats$Role,levels=c("Peripherals","Module hubs","Connectors","Network hubs"))
  pz <- ggplot(subset(stats,is.finite(Zi)),aes(Pi,Zi))+
    geom_vline(xintercept=.62,colour="grey65",linewidth=.3)+
    geom_hline(yintercept=2.5,colour="grey65",linewidth=.3)+
    geom_point(aes(colour=Role,shape=Role),size=1.15,alpha=.85)+
    scale_colour_manual(values=c(Peripherals="#B8B8B8",`Module hubs`="#26BEC4",Connectors="#FF1515",`Network hubs`="#762a83"),drop=TRUE)+
    scale_shape_manual(values=c(Peripherals=15,`Module hubs`=16,Connectors=17,`Network hubs`=18),drop=TRUE)+
    scale_x_continuous(breaks=seq(0,1,.2),expand=c(0,0))+
    scale_y_continuous(breaks=-2:3,expand=c(0,0))+
    coord_cartesian(xlim=c(0,1.1),ylim=range(c(-2.6,3.5,stats$Zi),na.rm=TRUE)*c(1.03,1.03))+
    annotate("text",x=c(.015,1.085,.015,1.085),y=c(3.32,3.32,-2.41,-2.41),
      label=c("Module hubs","Network connectors","Peripherals","Connectors"),hjust=c(0,1,0,1),size=3.1,family="serif")+
    labs(x=expression("Among module connectivities ("*P[i]*")"),
         y=expression("Within module connectivities ("*Z[i]*")"))+
    theme_bw(base_size=10,base_family="serif")+
    theme(panel.grid=element_blank(),panel.border=element_rect(colour="black",linewidth=.5),
      legend.position="inside",legend.position.inside=c(.84,.68),legend.title=element_blank(),
      legend.background=element_blank(),legend.key.size=unit(3.6,"mm"),legend.text=element_text(size=7),
      plot.margin=margin(2,2,2,2))
  role_plots[[group]] <- pz
  save_plot(pz,file.path(fig,paste0("13_",group,"_ZiPi_roles")),137,110)
}
write.csv(bind_rows(metrics),file.path(res,"13_visual_example_network_metrics.csv"),row.names=FALSE)

## 5. A/B/C/D 一张来源图保留一张组合图；不是拼接独立模板的 preview ====
shared_legend <- cowplot::get_legend(networks$T)
upper <- networks$G+theme(legend.position="none")+networks$T+theme(legend.position="none")+
  wrap_elements(shared_legend)+plot_layout(widths=c(1,1,.47))
lower <- role_plots$G+role_plots$T+plot_spacer()+plot_layout(widths=c(1,1,.47))
combined <- (upper/lower)+plot_layout(heights=c(1.15,.88))
combined <- combined+plot_annotation(tag_levels=list(c("A","B","","C","D",""))) &
  theme(plot.tag=element_text(size=13,face="bold",family="sans"),plot.tag.position=c(.01,.99))
save_plot(combined,file.path(fig,"13_networks_and_ZiPi"),265,202)
print(bind_rows(metrics))

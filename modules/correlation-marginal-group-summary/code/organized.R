# 相关性矩阵与边缘分布及分组摘要 / Correlation matrix with marginal bars and group summary
# 科学问题：代谢物之间的共变关系是否与其类别及平均丰度对应？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "correlation-marginal-group-summary")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：同步已确认的相关矩阵/边缘条形几何，保留本条目原始示例CSV与独立输出；输入/派生值均为合成示例，不冒充来源结果。
# 胆汁酸上三角相关矩阵、同坐标边缘浓度和五类汇总相关。
# 本轮依据 source/reference-02.png 恢复整图，而不是套用出版主题。
# 原始 CSV 是合成示例；下列结构化变换和来源标签不构成真实实验结果。
# 输入 data/bile_acid_*.csv；绘图输入及顺序另存 result/，不覆盖原 CSV。
Sys.setlocale("LC_CTYPE","English_United States.utf8")



source(file.path(module_root,"code/export.R"))
fig <- file.path(entry,"figure"); res <- file.path(entry,"result")
set.seed(260103)

## 1. 从现有合成观测构造可复用的结构化示例 ====
raw <- read.csv(file.path(module_root,"data/bile_acid_observations.csv"),row.names=1,check.names=FALSE)
# 来源图的标签顺序：不是声称这份合成数据的聚类顺序等同于论文。
labels <- c("Glycohyodeoxycholic acid","Arg-dihydroxylated BA","Arg-triHydroxylated BA",
  "Asn-dihydroxylated BA 1","Ser-triHydroxylated BA","Asn-triHydroxylated BA",
  "Ser-dihydroxylated BA","Cys-dihydroxylated BA","Ala-triHydroxylated BA","Cys-triHydroxylated BA",
  "Tauroursodeoxycholic acid","Tauro-\u03b1-Muricholic acid","Taurolithocholic acid",
  "Taurodeoxycholic acid","Taurocholic acid","Taurochenodeoxycholic acid",
  "Glycolithocholic acid","Glycodeoxycholic acid","Glycoursodeoxycholic acid",
  "Glycocholic acid","Glycochenodeoxycholic acid","Lithocholic acid","Deoxycholic acid",
  "Ursodeoxycholic acid","Cholic acid","Chenodeoxycholic acid","Phe-dihydroxylated BA",
  "Tyr-dihydroxylated BA","Tyrosocholic acid","Phenylalanocholic acid","Lys-dihydroxylated BA",
  "Lys-triHydroxylated BA","Met-triHydroxylated BA","Asn-dihydroxylated BA 2",
  "Gln-dihydroxylated BA","Trp-triHydroxylated BA","Leucocholic acid",
  "His-triHydroxylated BA","Glu-dihydroxylated BA","Gln-triHydroxylated BA","Glu-triHydroxylated BA")
group <- c(4,rep(1,9),4,2,4,4,2,2,4,4,4,2,2,5,5,5,3,3,rep(1,15))
n <- length(labels)
pcs <- scale(prcomp(scale(log10(raw+1e-6)))$x[,1:8])
e <- scale(log10(raw+1e-6))
visual <- matrix(0,nrow(raw),n,dimnames=list(rownames(raw),labels))
for(j in seq_len(n)) {
  # 共同因子恢复跨类别相关，分组因子保留模块；没有读取来源图的相关系数。
  module <- if(j<=10) 2 else if(j<=21) 3 else if(j<=26) 4 else 5
  common_weight <- if(j %in% c(22,23,24)) .15 else .55
  z <- common_weight*pcs[,1]+.58*pcs[,module]+.53*e[,(j-1)%%ncol(e)+1]
  target_logmean <- if(group[j]==1) 2.1+((j*7)%%13)/12 else 3.0+((j*11)%%17)/9
  visual[,j] <- 10^(target_logmean+.34*z)
}
r <- cor(log10(visual+1e-6))
group_cols <- c("#000000","#A69900","#8F1D00","#056258","#6A478F")
text_cols <- group_cols[group]
means <- log10(colMeans(visual))
sums <- sapply(1:5,function(g) rowSums(visual[,group==g,drop=FALSE]))
summary_r <- cor(log10(sums+1e-6))
write.csv(visual,file.path(res,"02_visual_example_abundance.csv"))
write.csv(data.frame(variable=labels,group=group,order=1:n,log10_mean=means),
          file.path(res,"02_visual_example_order.csv"),row.names=FALSE)
write.csv(r,file.path(res,"02_visual_example_correlations.csv"))
write.csv(summary_r,file.path(res,"02_visual_example_group_correlations.csv"))

## 2. 在同一毫米画布绘制矩阵、行标签与条形：不再让 patchwork 拉开基线 ====
pal <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
                         "#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(201)
cell <- 3.12; x0 <- 92; top <- 136; barx <- x0+n*cell+2
save_draw({
  par(mar=c(0,0,0,0),family="Arial",xaxs="i",yaxs="i",ps=10)
  plot.new(); plot.window(xlim=c(0,247),ylim=c(0,165),asp=1)
  for(i in 1:n) {
    y <- top-(i-.5)*cell
    for(j in i:n) {
      x <- x0+(j-.5)*cell
      rect(x-cell/2,y-cell/2,x+cell/2,y+cell/2,col="white",border="#B8B8B8",lwd=.45)
      symbols(x,y,circles=.46*cell*sqrt(abs(r[i,j])),inches=FALSE,
              bg=pal[round((r[i,j]+1)*100)+1],fg=NA,add=TRUE)
    }
    text(x0+(i-1)*cell-.9,y,labels[i],adj=c(1,.5),col=text_cols[i],cex=.62)
    text(x0+(i-.5)*cell,top+1.0,labels[i],adj=c(0,.5),srt=60,col=text_cols[i],cex=.62)
    # 对齐同一行中心；颜色表示绝对定量/微生物偶联，不是相关正负。
    rect(barx,y-.43*cell,barx+means[i]*3.6,y+.43*cell,
         col=if(group[i]==1) "#440154" else "#00A075",border=NA)
  }
  bottom <- top-n*cell
  segments(barx,bottom,barx+6*3.6,bottom,lwd=.65)
  segments(barx,bottom,barx,top,lwd=.65)
  for(tick in c(0,2,4,6)) {
    segments(barx+tick*3.6,bottom,barx+tick*3.6,bottom-1.3,lwd=.65)
    text(barx+tick*3.6,bottom-3.2,tick,cex=.78)
  }

  ## 3. 左侧四阶三角摘要及唯一的相关色标；恢复来源中的位置和完整圆点 ====
  sx <- 6; sy <- 115; sc <- 10.7
  for(i in 1:4) for(j in (i+1):5) {
    x <- sx+(j-1.5)*sc; y <- sy-(i-.5)*sc
    rect(x-sc/2,y-sc/2,x+sc/2,y+sc/2,col="white",border="#B8B8B8",lwd=.5)
    symbols(x,y,circles=.43*sc*sqrt(abs(summary_r[i,j])),inches=FALSE,
            bg=pal[round((summary_r[i,j]+1)*100)+1],fg=NA,add=TRUE)
  }
  for(j in 2:5) text(sx+(j-1.5)*sc,sy+3.4,j,cex=.9)
  for(i in 1:4) text(sx+(i-1)*sc-4,sy-(i-.5)*sc,i,cex=.9)
  for(k in 1:201) rect(6+(k-1)*43/201,sy-48,6+k*43/201,sy-44,col=pal[k],border=NA)
  rect(6,sy-48,49,sy-44,border="black",lwd=.6)
  for(tick in seq(-1,1,.2)) {
    x <- 6+(tick+1)*43/2
    segments(x,sy-48,x,sy-49,lwd=.5)
    text(x,sy-51,format(round(tick,1),nsmall=1),cex=.49)
  }
  legend_text <- c("1. Microbially conjugated bile acids","2. 1° conjugated bile acids",
                   "3. 1° bile acids","4. 2° conjugated bile acids","5. 2° bile acids")
  for(i in 1:5) {
    y <- 45-(i-1)*5.4
    rect(63,y-1.45,66,y+1.45,col=group_cols[i],border=NA)
    text(68,y,legend_text[i],adj=c(0,.5),col=group_cols[i],cex=.77)
  }
},file.path(fig,"02_correlation_marginals"),247,165)
cat("VISUAL_REPAIR: 41 variables; aligned upper-triangle grid and bars; synthetic inputs saved in result/\n")

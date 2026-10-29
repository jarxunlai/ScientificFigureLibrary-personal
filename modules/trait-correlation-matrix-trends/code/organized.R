# 水稻性状相关矩阵与年代趋势 / Rice trait matrix and year trends
# 科学问题：多个品质性状如何相互关联，各性状是否随育成年份变化？
# 输入：模块内合成CSV；只复用绘图层，不声称恢复原研究结论。
# 输出：显式传入目录的figure/和result/；未传目录时使用本次R会话临时目录。
# 运行：Rscript --vanilla code/organized.R <output-directory>
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
module_root <- dirname(dirname(normalizePath(script_file, winslash = "/")))
output_args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(output_args)) output_args[1] else file.path(tempdir(), "trait-correlation-matrix-trends")
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
root <- normalizePath(output_root, winslash = "/")
entry <- root
base <- root
draft <- root
batch <- root
for (output_dir in c("figure", "result")) dir.create(file.path(root, output_dir), showWarnings = FALSE, recursive = TRUE)
# 第二轮：收紧矩阵面板间隙、恢复右侧趋势的量纲与方向，统计仍从拟态输入计算；输入/派生值均为合成示例，不冒充来源结果。
Sys.setlocale("LC_CTYPE","English_United States.utf8")
# 微信相关性合集本地复刻；仅绘图层，输入为明确标记的合成示例。
# 输入：本条目 data/；输出：figure/*.png 和 *.pdf。




source(file.path(module_root,"code/export.R"))
fig <- file.path(entry, "figure")
dat <- file.path(module_root,"data")
theme_set(theme_classic(base_size = 9, base_family = "Arial"))
## 1. 读取水稻性状与年份 ====
d <- read.csv(file.path(dat, "rice_traits.csv"))
set.seed(270101)
z <- scale(d[,c("AAC","PC","LC","TV")])
d$AAC <- 19-.14*(d$Year-1970)+1.6*z[,1]
d$PC <- 8.4-.009*(d$Year-1970)+.28*z[,2]
d$LC <- 1.36-.0003*(d$Year-1970)+.055*z[,3]
d$TV <- 51+.32*(d$Year-1970)+4.2*z[,4]
write.csv(d,file.path(entry,"result/01_visual_example.csv"),row.names=FALSE)
vars <- c("AAC", "PC", "LC", "TV")
cols <- c(AAC = "#e6141c", PC = "#20a64f", LC = "#f4ba1a", TV = "#693290")
limits <- list(AAC=c(0,27), PC=c(0,12), LC=c(0,2.1), TV=c(0,90))
## 2. 矩阵与年代趋势；R平方与P值来自实际示例表 ====
save_base({
  layout(matrix(1:20, 4, 5, byrow=TRUE), widths=c(1,1,1,1,1.35))
  par(mar=c(.75,.75,.65,.2), mgp=c(1.1,.3,0), tcl=-.2, family="Arial", ps=11, cex=.9)
  for(i in 1:4) for(j in 1:5) {
    par(mar=c(if(i==4)2.2 else .5,if(j==5)2.8 else .65,.65,.3))
    if(j==5) {
      v <- vars[i]; fit <- lm(d[[v]] ~ d$Year)
      plot(d$Year,d[[v]],pch=16,cex=.7,col=cols[v],xlim=c(1970,2030),ylim=limits[[v]],xlab="",ylab="",xaxt=if(i==4) "s" else "n")
      abline(fit,lty=2,lwd=.7)
      text(1972, limits[[v]][2]*.14,sprintf("y = %.3f x %+.1f",coef(fit)[2],coef(fit)[1]),adj=0,cex=.8)
      text(1972, limits[[v]][2]*.045,sprintf("R^2=%.2f; P=%.2g",summary(fit)$r.squared,coef(summary(fit))[2,4]),adj=0,cex=.8)
      mtext(paste0(v," (%)"),side=2,line=1.3,cex=.75)
      if(i==4) mtext("Year",1,line=1.2,cex=.8)
    } else if(i==j) {
      v <- vars[i]; h <- hist(d[[v]],plot=FALSE,breaks=14)
      plot(h,freq=FALSE,col=cols[v],border="white",main="",axes=FALSE,xlab="",ylab="")
      curve(dnorm(x,mean(d[[v]]),sd(d[[v]])),add=TRUE,lwd=.8)
      box(col="grey55")
    } else if(i>j) {
      plot(d[[vars[j]]],d[[vars[i]]],pch=16,col="#bbbbbb99",cex=.7,xlab="",ylab="",axes=FALSE)
      abline(lm(d[[vars[i]]]~d[[vars[j]]]),col="#df7777",lty=2,lwd=.8); box(col="grey55")
    } else {
      ct <- cor.test(d[[vars[i]]],d[[vars[j]]]); pv <- ct$p.value
      bg <- if(pv<.001) "#FFDDDD" else if(pv<.01) "#FFEEEE" else if(pv<.05) "#FFFFCC" else "#EEEEEE"
      plot.new(); rect(0,0,1,1,col=bg,border="grey55")
      text(.5,.5,sprintf("R² = %.3f\nP = %.2g\n%s",unname(ct$estimate)^2,pv,if(pv<.001) "***" else if(pv<.01) "**" else if(pv<.05) "*" else "ns"),cex=1.0)
    }
    if(j<=4 && i==4) mtext(vars[j],1,line=.6,cex=.85)
    if(j==1) mtext(vars[i],2,line=.6,cex=.85)
  }
},file.path(fig,"01_rice_matrix_year_trends"),220,145)

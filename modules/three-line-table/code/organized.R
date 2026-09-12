# =============================================================================
# 多组哑铃图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# =============================================================================

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
set.seed(20260912)
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(grid)
# 一、模拟观测，不包含真实医院或个体记录。
linelist <- data.frame(hospital=rep(paste("Hospital",LETTERS[1:6]),each=120),outcome=sample(c("Death","Recover"),720,TRUE),ct_blood=round(rnorm(720,22,4),1))
# 二、各组和总计使用同一纳入范围。
summary_rows <- linelist %>% group_by(hospital,outcome) %>% summarise(N=n(),ct=median(ct_blood),.groups="drop")
total_rows <- linelist %>% group_by(outcome) %>% summarise(N=n(),ct=median(ct_blood),.groups="drop") %>% mutate(hospital="Total")
tab <- bind_rows(summary_rows,total_rows) %>% pivot_wider(names_from=outcome,values_from=c(N,ct)) %>% mutate(Known=N_Death+N_Recover,Pct_Recover=sprintf("%.1f%%",100*N_Recover/Known),Pct_Death=sprintf("%.1f%%",100*N_Death/Known)) %>% select(hospital,Known,N_Recover,Pct_Recover,ct_Recover,N_Death,Pct_Death,ct_Death)
stopifnot(tail(tab$Known,1)==sum(head(tab$Known,-1)))
# 三、三条主水平线、双层表头与总计加粗。
png(file.path(root,"preview.png"),width=2000,height=900,res=180,type="cairo")
grid.newpage()
xs <- c(.1,.26,.38,.49,.60,.71,.82,.93)
grid.text("Illustrative summary (synthetic data)",.5,.96,gp=gpar(fontsize=14,fontface="bold"))
grid.segments(.02,.90,.98,.90,gp=gpar(lwd=2))
grid.text("Hospital",xs[1],.82,gp=gpar(fontface="bold"))
grid.text("Known\noutcomes",xs[2],.82,gp=gpar(fontface="bold"))
grid.text("Recovered",mean(xs[3:5]),.85,gp=gpar(fontface="bold"))
grid.text("Died",mean(xs[6:8]),.85,gp=gpar(fontface="bold"))
grid.text(rep(c("Total","Percent","Median CT"),2),xs[3:8],.77,gp=gpar(fontsize=10))
grid.segments(.02,.71,.98,.71,gp=gpar(lwd=1.2))
for(i in seq_len(nrow(tab))) {
 y <- .66-(i-1)*.078
 for(j in seq_len(ncol(tab))) grid.text(as.character(tab[i,j][[1]]),xs[j],y,gp=gpar(fontsize=11,fontface=ifelse(tab$hospital[i]=="Total","bold","plain")))
}
grid.segments(.02,.14,.98,.14,gp=gpar(lwd=2))
dev.off()

# =============================================================================
# 批量散点拟合图
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
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})

# 数据构建：
data <- data.frame(MTAP = sort(runif(50, min = 4, max = 6))[c(1:8,sample(9:41),42:50)],
           PRMT1 = sort(runif(50, min = 11, max = 14), decreasing = T)[c(1:8,sample(9:41),42:50)])

# library(Hmisc)  # 当前 Pixi R 下该包原生 DLL 会崩溃；使用 shim rcorr

# 计算相关系数和p值
res <- rcorr(data$MTAP, data$PRMT1)
p_value <- res$P[1,2]
cor_value <- round(res$r[1,2], 2)

# 绘图：
ggplot(data,aes(MTAP, PRMT1))+
  geom_point(color = "#988d7b")+
  geom_smooth(method = "lm", formula = y ~ x, 
              # 调整置信区间颜色：
              fill = "#b2e7fa", color = "#00aeef", alpha = 0.8)+
  theme_bw()+
  theme(
    # 去除网格线：
    panel.grid = element_blank(),
    # 修改坐标轴标签
    axis.title = element_text(face = "bold.italic"),
    # 标题居中：
    plot.title = element_text(hjust = 0.5)
    )+
  labs(title = paste0("ρ =", cor_value, ", q = ", p_value))


# 数据新增5列：
data$NECTN2 <- sort(runif(50, min = 13, max = 16), decreasing = T)[c(1:8,sample(9:41),42:50)]
data$IDO1 <- sort(runif(50, min = 12, max = 18), decreasing = T)[c(1:8,sample(9:41),42:50)]
data$SIRPA <- sort(runif(50, min = 10, max = 16), decreasing = T)[c(1:8,sample(9:41),42:50)]
data$SIRPA2<- sort(runif(50, min = 9, max = 14), decreasing = T)[c(1:8,sample(9:41),42:50)]
data$MIF <- sort(runif(50, min = 11, max = 17), decreasing = T)[c(1:8,sample(9:41),42:50)]

# 循环计算相关系数并绘图：
# 创建一个空列表，存储返回的图形：
p_list <- list()

for (i in 2:ncol(data)) {
  res <- rcorr(data$MTAP, data[,i])
  p_value <- signif(res$P[1,2], 2)
  cor_value <- round(res$r[1,2], 2)
  
  # 每次新建一个绘图数据框：
  data_new <- data[,c(1,i)]
  colnames(data_new) <- c("MTAP", "y")
  
  p <- ggplot(data_new,aes(x = MTAP, y = y))+
    geom_point(color = "#988d7b")+
    geom_smooth(method = "lm", formula = y ~ x, 
                # 调整置信区间颜色：
                fill = "#b2e7fa", color = "#00aeef", alpha = 0.8)+
    theme_bw()+
    ylab(colnames(data)[i])+
    theme(
      # 去除网格线：
      panel.grid = element_blank(),
      # 修改坐标轴标签
      axis.title = element_text(face = "bold.italic"),
      # 标题居中：
      plot.title = element_text(hjust = 0.5, size = 10)
    )+
    labs(title = paste0("r =", cor_value, ", q = ", p_value))
  p_list[[i-1]] <- p
}

# 拼图
library(cowplot)
library(patchwork)

p <- plot_grid(p_list[[1]], p_list[[2]], p_list[[3]],
          p_list[[4]], p_list[[5]], p_list[[6]], ncol = 3)

ggsave(file.path(out_dir, "plot.pdf"), plot = p, height = 6, width = 9)

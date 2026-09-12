# =============================================================================
# 复杂网络图
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

library(readxl)
# library(tidyverse) 已在文件头拆包
library(ggplot2)

# 最小关系表，按实际节点数构建布局，不再硬编码 120 个不同节点。
data <- read.csv(file.path(root,"data","input.csv"),check.names=FALSE)


# 统计代谢物及其分类出现的次数，以确定散点的大小：
metabolite_count <- as.data.frame(table(data$`Serum metabolite`))
Polutant_count <- as.data.frame(table(data$Pollutant))
Outcome_count <- as.data.frame(table(data$Outcome))

# 给每个散点赋予一个坐标轴位置：
metabolite_count$x[1:nrow(metabolite_count)] <- 1:nrow(metabolite_count)
metabolite_count$y[1:nrow(metabolite_count)] <- 0
# 这个数字根据后面的图形自行修改：
Polutant_count$x[1:nrow(Polutant_count)] <- seq(1,nrow(metabolite_count),length.out=nrow(Polutant_count))
Polutant_count$y[1:nrow(Polutant_count)] <- 5
# 这个数字根据后面的图形自行修改：
Outcome_count$x[1:nrow(Outcome_count)] <- seq(1,nrow(metabolite_count),length.out=nrow(Outcome_count))
Outcome_count$y[1:nrow(Outcome_count)] <- -5

# 合并数据：
data_count <- rbind(rbind(metabolite_count, Polutant_count), Outcome_count)
colnames(data_count)[1] <- "Serum metabolite"
# 加入分组信息 -- 到这绘制散点的数据算是完成了！
data_count <- left_join(data_count, unique(data[,1:2]), by = "Serum metabolite")

idx_m <- seq_len(nrow(metabolite_count))
idx_p <- nrow(metabolite_count)+seq_len(nrow(Polutant_count))
idx_o <- nrow(metabolite_count)+nrow(Polutant_count)+seq_len(nrow(Outcome_count))
data_count$`Super pathway`[idx_p] <- "OPEs"
data_count$`Super pathway`[idx_o] <- data_count$`Serum metabolite`[idx_o]

########## 先绘制散点 ----------------
colors <- c("#476b71", "#8697a0", "#2da3d1", "#806766",
            "#55bfe2", "#b2bec5", "#d2b698")
names(colors) <- unique(data_count$`Super pathway`)[1:7]
colors <- c(colors, "OPEs" = "#1999a9", "GSP" = "#efc000",
            "HOMA-IR" = "#9a8419", "FPG" = "#1981c8")

p <- ggplot(data_count)+
  geom_point(aes(x, y, size = Freq, color = `Super pathway`))+
  geom_text(data = data_count[idx_m,],
            aes(x, y-0.1, label = `Serum metabolite`, color = `Super pathway`),
            angle = 90, hjust = 1, vjust = 0.5, size = 1.5, show.legend = F)+
  geom_text(data = data_count[idx_p,],
            aes(x+3, y, label = `Serum metabolite`, color = `Super pathway`),
            angle = 0, hjust = 0, size = 4, show.legend = F)+
  geom_text(data = data_count[idx_o,],
            aes(x, y-0.5, label = `Serum metabolite`, color = `Super pathway`),
            angle = 0, hjust = 0.5, vjust = 0.5, size = 4, show.legend = F)+
  scale_color_manual(name = "Class", values = colors)+
  theme_void()

p

########### 折线的数据 ---------------
data_line <- data[, c(1,3)]
data_line2 <- data[,c(1,4)]
colnames(data_line2) <- colnames(data_line)
data_line <- rbind(data_line,data_line2)
data_line$group <- paste0("group", 1:nrow(data_line))
data_line <- pivot_longer(data_line, cols = -group,
                          names_to = "Class", values_to = "Serum metabolite")

data_line <- left_join(data_line, unique(data_count[,c(1,3,4)]), by = "Serum metabolite")

final_plot <- p+geom_line(data = data_line, aes(x, y, group = group), colour = "#b7bfcb",
            linewidth = 0.2, alpha = 0.3, show.legend = F)

ggsave(file.path(out_dir, "plot.pdf"), height = 5, width = 10)

ggplot2::ggsave(file.path(root,"preview.png"),plot=final_plot,width=12,height=6,dpi=200,bg="white")

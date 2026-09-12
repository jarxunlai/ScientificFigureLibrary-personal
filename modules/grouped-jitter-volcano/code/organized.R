# =============================================================================
# 单细胞分组抖动火山图
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
# 数据路径：file.path(root, "data", ...)
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 读取数据：
data <- read.csv(file.path(root, "data", "data02.csv"), row.names = 1)

# 选取上下调超过25%，且矫正p值小于0.05：
padj <- 0.05
data$change <- ifelse(data$log2FoldChange >= log2(1.25) & data$padj < padj, 
                      "Up regulate",
                      ifelse(data$log2FoldChange < log2(0.75) & data$padj < padj,
                             "Down regulate", 
                             "Not significant"))

# 创建空列表，用于存放每组的差异分析数据：
# 这种创建方法只是创建个数据，自己用的时候没有参考价值，根据自己的数据情况调整：
DEG_list <- list()
for (i in 1:8) {
  data_tmp <- data[sample(1:nrow(data), 2000), ]
  data_tmp$group <- paste("cell", LETTERS[i], sep = "_")
  DEG_list[[paste0("cell", LETTERS[i])]] <- data_tmp
}

plot_dat <- dplyr::bind_rows(DEG_list)

# 这里只选取log2FC绝对值排前20的gene；
label_data <- plot_dat[plot_dat$change != "Not significant", ]
label_data <- label_data[order(label_data$log2FoldChange, decreasing = T)[c(1:10, (nrow(label_data)-10):nrow(label_data))], ]

label_data$label <- rownames(label_data)

# 美化1：
ggplot(plot_dat)+
  # 抖动散点：
  geom_jitter(aes(group, log2FoldChange, color = change),
              size=0.85, width = 0.4, alpha= .8)+
  # 分组方块：
  geom_tile(aes(group, 0, fill = group),
            height=0.8,
            color = "black",
            alpha = 0.5,
            show.legend = F,
            width=0.85) +
  # 文字：
  geom_text(data = plot_dat[!duplicated(plot_dat$group), ], 
            aes(group, 0, label = group),
            size =2,
            color ="black") +
  # 基因标签：
  geom_text_repel(
    data = label_data,
    aes(group, log2FoldChange, label = label),
    size=2, max.overlaps = 100
  ) +
  xlab("Cell Subtype")+
  ylab("log2FoldChange")+
  # 颜色模式
  scale_fill_manual(values = brewer.pal(8, "Set3"))+
  scale_color_manual(name = "Regulate", values = c("#b3de69","#999999","#fb8072"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "top")

ragg::agg_png(file.path(out_dir, "grouped_volcano_p1.png"), width = 7, height = 7, units = "in", res = 300)
print(last_plot())
dev.off()
file.copy(file.path(out_dir, "grouped_volcano_p1.png"), file.path(root, "preview.png"), overwrite = TRUE)
ggsave(file.path(out_dir, "grouped_volcano_p1.pdf"), last_plot(), width = 7, height = 7)


# 美化2：
ggplot(plot_dat)+
  # 彩色散点：
  geom_jitter(data = plot_dat[plot_dat$change != "Not significant",], 
              aes(group, log2FoldChange, color = group),
              size=0.85, width = 0.4, alpha= .8)+
  # 灰色散点
  geom_jitter(data = plot_dat[plot_dat$change == "Not significant",], 
              aes(group, log2FoldChange), color = "#999999",
              size=0.85, width = 0.4, alpha= .8)+
  # 分组方块：
  geom_tile(aes(group, 0, fill = group),
            height=0.8,
            color = "black",
            alpha = 0.5,
            show.legend = F,
            width=0.85) +
  # 文字：
  geom_text(data = plot_dat[!duplicated(plot_dat$group), ], 
            aes(group, 0, label = group),
            size =2,
            color ="black") +
  # 基因标签：
  geom_text_repel(
    data = label_data,
    aes(group, log2FoldChange, label = label),
    size=2, max.overlaps = 100
  )+
  xlab("Cell Subtype")+
  ylab("log2FoldChange")+
  # 颜色模式
  scale_fill_manual(values = brewer.pal(8, "Set3"))+
  scale_color_manual(name = "Cell subtype", values = brewer.pal(8, "Set3"))+
  #scale_color_manual(name = "Regulate", values = c("#b3de69","#999999","#fb8072"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "top")

ragg::agg_png(file.path(out_dir, "grouped_volcano_p2.png"), width = 7, height = 7, units = "in", res = 300)
print(last_plot())
dev.off()
ggsave(file.path(out_dir, "grouped_volcano_p2.pdf"), last_plot(), width = 7, height = 7)


  
  
  
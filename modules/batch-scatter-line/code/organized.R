# =============================================================================
# 批量散点与折线图
# organized：线性脚本 + 中文分节；路径相对本条目目录
# =============================================================================
# 来源：KS科研分享「生信绘图」042批量散点图+批量折线图
# 本地复现：drafts/ks-shengxin-huitu-repro/042-batch-scatter-line
# Pixi：项目根 default 环境；library(tidyverse) 已拆成 ggplot2/dplyr/tidyr 等。
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

# 载入R包：
library(ggplot2)
# library(tidyverse) 已在文件头拆包

# 批量构建数据：
for (i in 1:6) {
  base_num <- sample(2:6, 1)
  data <- data.frame(x = rep(0:19, 2),
                     y = c(log(1:20, base = base_num)*1000,
                           log(1:20, base = base_num)*900),
                     group = rep(c("group1", "group2"), each = 20)
  )
  assign(paste0("data", i), data)
}

############ 批量绘制散点图 ------------
fills <- c("#e5f5f0", "#e5f0f7", "#edf7fc", "#fefdf4", "#fcf5e5", "#f9f1f6")
colors <- c("#0c7d5e", "#0c5e8c", "#56b4e9", "#918912", "#b68518", "#b14481")

p_list <- list()
anno_texts <- c("Half precip.\n and no clipping",
                "Normal precip.\n and no clipping",
                "Double precip.\n and no clipping",
                "Half precip.\n and clipping",
                "Normal precip.\n and clipping",
                "Double precip.\n and clipping")

for (i in 1:6) {
  if (i == 1) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_point(aes(x, y, color = group), size = 1)+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            axis.text.x = element_blank())
  } else if (i %in%  c(2:3)) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_point(aes(x, y, color = group), size = 1)+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
  } else if (i == 4) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_point(aes(x, y, color = group), size = 1)+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white")
            # axis.text.x = element_blank()
            )
  } else {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_point(aes(x, y, color = group), size = 1)+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            # axis.text.x = element_blank()
            axis.text.y = element_blank()
            )
  }
}

library(cowplot)

p1 <- plot_grid(plotlist = p_list, ncol = 3, rel_widths = c(1.1, 1, 1))

ggsave(file.path(out_dir, "scatter_plot.pdf"), plot = p1, height = 4, width = 8)


############ 批量绘制散点图 ------------
fills <- c("#e5f5f0", "#e5f0f7", "#edf7fc", "#fefdf4", "#fcf5e5", "#f9f1f6")
colors <- c("#0c7d5e", "#0c5e8c", "#56b4e9", "#918912", "#b68518", "#b14481")

p_list <- list()
anno_texts <- c("Half precip.\n and no clipping",
                "Normal precip.\n and no clipping",
                "Double precip.\n and no clipping",
                "Half precip.\n and clipping",
                "Normal precip.\n and clipping",
                "Double precip.\n and clipping")

for (i in 1:6) {
  if (i == 1) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_line(aes(x, y, color = group))+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            axis.text.x = element_blank())
  } else if (i %in%  c(2:3)) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_line(aes(x, y, color = group))+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank())
  } else if (i == 4) {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_line(aes(x, y, color = group))+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white")
            # axis.text.x = element_blank()
      )
  } else {
    p_list[[i]] <- ggplot(get(paste0("data", i)))+
      geom_line(aes(x, y, color = group))+
      scale_color_manual(values = c("#0072b2", "#d55e00"))+
      xlab("")+
      ylab("")+
      annotate("text", x = 20, y = 0, label = anno_texts[i],
               color = colors[i], hjust = 1, vjust = 0)+
      theme_bw()+
      theme(legend.position = "none",
            panel.background = element_rect(fill = fills[i]),
            panel.grid = element_line(color = "white"),
            # axis.text.x = element_blank()
            axis.text.y = element_blank()
      )
  }
}

library(cowplot)

p2 <- plot_grid(plotlist = p_list, ncol = 3, rel_widths = c(1.1, 1, 1))

ggsave(file.path(out_dir, "line_plot.pdf"), plot = p2, height = 4, width = 8)


#### 再拼图：
p3 <- plot_grid(p2, p1, labels = c("e", "f"))

ggsave(file.path(out_dir, "all_plot.pdf"), plot = p3, height = 4, width = 14)

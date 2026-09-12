# =============================================================================
# 批量火山图
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
library(cowplot)
library(grid)

# 模拟数据：仅演示多组比较布局，不对应真实效应或统计结论。
Pollutant_names <- c("TiBP", "TMPP", "TPHP", "DBP", "DPHP")
data <- data.frame(Pollutant = factor(rep(Pollutant_names, each = 200), levels = Pollutant_names),
 Estimate = rnorm(1000, 0, 0.12), FDR = 10^(-runif(1000, 0, 5)), Assiciation = "Not significant")

# 新建group列：如果FDR不显著，为灰色；
data$group <- data$Assiciation
data$group[data$FDR<0.05 & data$Estimate<0] <- "Negative"
data$group[data$FDR<0.05 & data$Estimate>0] <- "Positive"
data$group[data$FDR>0.05] <- "Not significant"

# 单独绘图：
plot_data <- data %>% filter(Pollutant == Pollutant_names[1])
xlim_value <- max(abs(plot_data$Estimate*100))
ylim_value <- max(-log10(plot_data$FDR))
p1 <- ggplot(plot_data)+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = -log10(0.05), linetype = 3)+
  geom_point(aes(Estimate*100, -log10(FDR), size = abs(Estimate),
                 fill = group), shape = 21, color = "white")+
  annotate("text", x = -xlim_value*3/4, y = ylim_value*7/8, color = "#50b0d4",
           label = paste0("Negative:", sum(plot_data$group == "Negative")))+
  annotate("text", x = xlim_value*3/4, y = ylim_value*7/8, color = "#d14d49",
           label = paste0("Positive:", sum(plot_data$group == "Positive")))+
  xlab("Simulated effect (scaled)")+
  scale_fill_manual(values = c("Negative" = "#50b0d4",
                               "Positive" = "#d14d49",
                               "Not significant" = "#a8b1ae"))+
  scale_size_continuous(range = c(1, 3))+
  xlim(-xlim_value, xlim_value)+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# 由于标题有个背景，ggplot2不太好加，分面的话又会导致坐标轴不好加，
# 这个办法也是权宜之计，属实是有点杀鸡用牛刀了，可以在AI中轻松绘制；
title_plot <- ggplot() +
  theme_void() +
  xlim(-0.108, 0.092)+
  annotation_custom(grob = rectGrob(gp = gpar(fill = "#2495a4", col = NA)),
                    xmin = -0.1, xmax = 0.1, ymin = -0.5, ymax = 1.5)+
  annotate("text", x = 0, y = 0.5, label = Pollutant_names[1],
           size = 5, hjust = 0.5, color = "white")

pdf(file.path(out_dir, "plot.pdf"), height = 4, width = 5)
plot_grid(title_plot, p1, ncol = 1, rel_heights = c(1, 10))
dev.off()



############ 批量绘图 -------------
plot_list <- list()
for (i in 1:length(Pollutant_names)) {
  plot_data <- data %>% filter(Pollutant == Pollutant_names[i])
  xlim_value <- max(abs(plot_data$Estimate*100))
  ylim_value <- max(-log10(plot_data$FDR))
  if (i == 1) {
    p1 <- ggplot(plot_data)+
      geom_vline(xintercept = 0, linetype = 3)+
      geom_hline(yintercept = -log10(0.05), linetype = 3)+
      geom_point(aes(Estimate*100, -log10(FDR), size = abs(Estimate),
                     fill = group), shape = 21, color = "white")+
      annotate("text", x = -xlim_value*2/3, y = ylim_value*7/8, color = "#50b0d4",
               label = paste0("Negative:", sum(plot_data$group == "Negative")))+
      annotate("text", x = xlim_value*2/3, y = ylim_value*7/8, color = "#d14d49",
               label = paste0("Positive:", sum(plot_data$group == "Positive")))+
      xlab("")+
      scale_fill_manual(values = c("Negative" = "#50b0d4",
                                   "Positive" = "#d14d49",
                                   "Not significant" = "#a8b1ae"))+
      scale_size_continuous(range = c(1, 3))+
      xlim(-xlim_value, xlim_value)+
      theme_classic()+
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))


    # 由于标题有个背景，ggplot2不太好加，分面的话又会导致坐标轴不好加，
    # 这个办法也是权宜之计，属实是有点杀鸡用牛刀了，可以在AI中轻松绘制；
    title_plot <- ggplot() +
      theme_void() +
      xlim(-0.12, 0.08)+
      annotation_custom(grob = rectGrob(gp = gpar(fill = "#2495a4", col = NA)),
                        xmin = -0.1, xmax = 0.1, ymin = -0.5, ymax = 1.5)+
      annotate("text", x = 0, y = 0.5, label = Pollutant_names[i],
               size = 5, hjust = 0.5, color = "white")
  } else {
    p1 <- ggplot(plot_data)+
      geom_vline(xintercept = 0, linetype = 3)+
      geom_hline(yintercept = -log10(0.05), linetype = 3)+
      geom_point(aes(Estimate*100, -log10(FDR), size = abs(Estimate),
                     fill = group), shape = 21, color = "white")+
      annotate("text", x = -xlim_value*2/3, y = ylim_value*7/8, color = "#50b0d4",
               label = paste0("Negative:", sum(plot_data$group == "Negative")))+
      annotate("text", x = xlim_value*2/3, y = ylim_value*7/8, color = "#d14d49",
               label = paste0("Positive:", sum(plot_data$group == "Positive")))+
      xlab("")+
      ylab("")+
      scale_fill_manual(values = c("Negative" = "#50b0d4",
                                   "Positive" = "#d14d49",
                                   "Not significant" = "#a8b1ae"))+
      scale_size_continuous(range = c(1, 3))+
      xlim(-xlim_value, xlim_value)+
      theme_classic()+
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))


    # 由于标题有个背景，ggplot2不太好加，分面的话又会导致坐标轴不好加，
    # 这个办法也是权宜之计，属实是有点杀鸡用牛刀了，可以在AI中轻松绘制；
    title_plot <- ggplot() +
      theme_void() +
      xlim(-0.12, 0.08)+
      annotation_custom(grob = rectGrob(gp = gpar(fill = "#2495a4", col = NA)),
                        xmin = -0.1, xmax = 0.1, ymin = -0.5, ymax = 1.5)+
      annotate("text", x = 0, y = 0.5, label = Pollutant_names[i],
               size = 5, hjust = 0.5, color = "white")
  }

  plot_list[[i]] <- plot_grid(title_plot, p1, ncol = 1, rel_heights = c(1, 10))
}

pdf(file.path(out_dir, "plot_all.pdf"), height = 3, width = 15)
plot_grid(plotlist = plot_list, nrow = 1, rel_widths = c(1.1, 1, 1, 1, 1))
dev.off()

ggplot2::ggsave(file.path(root, "preview.png"), plot = cowplot::plot_grid(plotlist = plot_list, nrow = 1, rel_widths = c(1.1,1,1,1,1)), width = 15, height = 3, dpi = 200, bg = "white")

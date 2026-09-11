# =============================================================================
# 柱状抖动散点误差棒
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

# library(tidyverse) 已在文件头拆包

# 构造模拟数据：
for (i in 2010:2016) {
  data <- data.frame(A = sample(-1500:1000, 5),
                     B = sample(-1000:100, 5),
                     C = sample(-1500:10, 5),
                     D = sample(-1000:0, 5),
                     E = sample(-1000:100, 5),
                     F = sample(-1000:100, 5))

  # 宽数据转长数据：
  data_long <- pivot_longer(data, cols = everything(),
                            names_to = "group", values_to = "Richness") %>%
    group_by(group) %>%
    mutate(means = mean(Richness))

  assign(paste0("data_", i, "_long"), data_long)
}


# 绘制单图：
ggplot(data_2010_long)+
  # errorbar：
  stat_boxplot(aes(group, Richness, color = group), geom = "errorbar", width=0.2)+
  # 柱状图：
  geom_col(data=unique(data_2010_long[,c(1,3)]),
           aes(group, means, fill = group, color = group), alpha = 0.6)+
  # 抖动散点图：
  geom_jitter(aes(group, Richness, color = group), width = 0.2)+
  # 颜色模式：
  scale_fill_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                               "#d5c711", "#af7a06", "#b24581"))+
  scale_color_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                               "#d5c711", "#af7a06", "#b24581"))+
  ylim(-2000, 1500)+
  ggtitle("2010")+
  # 主题
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_line(color = "white"),
        panel.background = element_rect(fill = "#fff6e3"))

ggsave(file.path(out_dir, "single_plot.pdf"), height = 4, width = 4.5)


# 循环作图：
# 背景色由浅变深：
fills <- colorRampPalette(c("#fffcf5", "#fff0cd"))(7)
p_list <- list()

m = 1
for (i in 2010:2016) {
  if (i == 2010|i == 2015) {
    p <- ggplot(get(paste0("data_", i, "_long")))+
      # errorbar：
      stat_boxplot(aes(group, Richness, color = group), geom = "errorbar", width=0.2)+
      # 柱状图：
      geom_col(data=unique(get(paste0("data_", i, "_long"))[,c(1,3)]),
               aes(group, means, fill = group, color = group), alpha = 0.6)+
      # 抖动散点图：
      geom_jitter(aes(group, Richness, color = group), width = 0.2)+
      # 颜色模式：
      scale_fill_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                                   "#d5c711", "#af7a06", "#b24581"))+
      scale_color_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                                    "#d5c711", "#af7a06", "#b24581"))+
      ylim(-2000, 1500)+
      xlab("")+
      ylab("")+
      ggtitle(i)+
      # 主题
      theme_bw()+
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_line(color = "white"),
            panel.background = element_rect(fill = fills[m]))
    p_list[[m]] <- p
    assign(paste0("p_",m), p)
    m = m+1
  } else {
    p <- ggplot(get(paste0("data_", i, "_long")))+
      # errorbar：
      stat_boxplot(aes(group, Richness, color = group), geom = "errorbar", width=0.2)+
      # 柱状图：
      geom_col(data=unique(get(paste0("data_", i, "_long"))[,c(1,3)]),
               aes(group, means, fill = group, color = group), alpha = 0.6)+
      # 抖动散点图：
      geom_jitter(aes(group, Richness, color = group), width = 0.2)+
      # 颜色模式：
      scale_fill_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                                   "#d5c711", "#af7a06", "#b24581"))+
      scale_color_manual(values = c("#31ae88", "#0072b2", "#2997d6",
                                    "#d5c711", "#af7a06", "#b24581"))+
      xlab("")+
      ylab("")+
      ylim(-2000, 1500)+
      ggtitle(i)+
      # 主题
      theme_bw()+
      theme(legend.position = "none",
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_line(color = "white"),
            panel.background = element_rect(fill = fills[m]))
    p_list[[m]] <- p
    assign(paste0("p_",m), p)
    m = m+1
  }
}

# 拼图
library(cowplot)

plot_grid(plotlist = p_list, ncol = 5, rel_widths = c(1.1, 1, 1, 1, 1))

ggsave(file.path(out_dir, "plot_all.pdf"), height = 6, width = 15)

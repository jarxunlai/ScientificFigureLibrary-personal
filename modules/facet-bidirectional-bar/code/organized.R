# =============================================================================
# 分面双向柱加热图
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

# library(tidyverse) 已在文件头拆包
library(patchwork)

# 构造模拟数据：
Virome <- runif(61, 0, 0.2)
Bacteriome <- runif(61, 0, 1)

group <- factor(rep(c("Anthropometric", "Lifestyle", "Diet", "Disease",
             "BSS", "PhysicalActivity", "Medication"),
             c(5, 2, 16, 12, 5, 1, 20)),
             levels = c("Anthropometric", "Lifestyle", "Diet", "Disease",
                        "BSS", "PhysicalActivity", "Medication"))

bar_data <- data.frame(Virome = Virome,
                       Bacteriome = Bacteriome,
                       group = group)

# 排序：
bar_data_sort1 <- bar_data %>%
  group_by(group) %>%
  mutate(Virome = sort(Virome, decreasing = T))

bar_data_sort2 <- bar_data %>%
  group_by(group) %>%
  mutate(Bacteriome = sort(Bacteriome, decreasing = T))

p1 <- ggplot(bar_data_sort1)+
  annotate("rect", xmin = seq(1.5, 59.5, 2),
           xmax = seq(2.5, 60.5, 2), ymin = 0, ymax = 0.2,
           fill = rep("#f5f5f5", 30))+
  geom_hline(yintercept = seq(0.05, 0.2, 0.05), color = "#e6e6e6")+
  geom_col(aes(x = 1:nrow(bar_data_sort1), y = Virome, fill = group))+
  scale_fill_manual(name = "",values = c("#4f6980", "#849db1", "#a2ceaa", "#638b66",
                      "#bfbb60", "#f47942", "#fbb04e"))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  xlab("")+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        panel.grid = element_blank())+
  guides(fill = guide_legend(nrow = 1))

p1

ggsave(file.path(out_dir, "p1.pdf"), height = 3, width = 12)

p2 <- ggplot(bar_data_sort2)+
  annotate("rect", xmin = seq(1.5, 59.5, 2),
           xmax = seq(2.5, 60.5, 2), ymin = 0, ymax = 1,
           fill = rep("#f5f5f5", 30))+
  geom_hline(yintercept = seq(0.25, 1, 0.25), color = "#e6e6e6")+
  geom_col(aes(x = 1:nrow(bar_data_sort1), y = Bacteriome, fill = group))+
  scale_fill_manual(values = c("#4f6980", "#849db1", "#a2ceaa", "#638b66",
                               "#bfbb60", "#f47942", "#fbb04e"))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  xlab("")+
  theme_bw()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+
  scale_y_reverse()

p2

ggsave(file.path(out_dir, "p2.pdf"), height = 3, width = 12)

p1/p2

ggsave(file.path(out_dir, "p1_p2.pdf"), height = 6, width = 12)


# 矩阵热图数据：
p_mat <- matrix(NA, nrow = 2, ncol = 61)
colnames(p_mat) <- paste0("feature_", 1:61)

anno_mat <- matrix(NA, nrow = 2, ncol = 61)
sig_index <- sample(1:122, 20)
p_mat[sig_index] <- runif(20, -20, 20)
anno_mat[sig_index] <- "*"

library(ComplexHeatmap)
library(circlize)

pdf(file.path(out_dir, "p3.pdf"), height = 1.5, width = 12)
Heatmap(p_mat,
        col = colorRamp2(c(-20, 0, 20),
                         c("#69a9d2", "white", "#c44c4b")),
        na_col = "white",
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(anno_mat[i, j])){
            grid.text(sprintf("*"),
                      x, y, gp = gpar(fontsize = 10))}
        },
        rect_gp = gpar(col = "black"),
        cluster_rows = F,
        cluster_columns = F,
        show_heatmap_legend = F)

dev.off()

ggplot2::ggsave(file.path(root, "preview.png"), plot = p1/p2, width = 12, height = 6, dpi = 200, bg = "white")

# 技术性导出工具；布局、颜色、拟合和阈值均留在条目脚本。
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)
options(stringsAsFactors = FALSE)
options(device = function(...) ragg::agg_capture(width=7,height=7,units="in",res=120,background="white"))
save_plot <- function(plot, stem, width = 120, height = 100, dpi = 300) {
  ggsave(paste0(stem, ".png"), plot = plot, device = ragg::agg_png,
         width = width, height = height, units = "mm", dpi = dpi, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot = plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "mm", family = "Arial", bg = "white")
}

# 仅负责设备和路径；统计选择、布局和颜色保留在各绘图脚本。
options(stringsAsFactors = FALSE)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(grid)
options(device = function(...) ragg::agg_capture(width=7,height=7,units="in",res=120,background="white"))
save_plot <- function(plot, stem, width = 180, height = 150, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, device = ragg::agg_png,
         width = width, height = height, units = "mm", dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "mm", bg = "white", family = "Arial")
}
save_base <- function(expr, stem, width = 180, height = 150, dpi = 300) {
  e <- substitute(expr); env <- parent.frame()
  ragg::agg_png(paste0(stem, ".png"), width = width, height = height, units = "mm", res = dpi, background = "white")
  eval(e, envir = env); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width/25.4, height = height/25.4, family = "Arial", bg = "white")
  eval(e, envir = env); dev.off()
}

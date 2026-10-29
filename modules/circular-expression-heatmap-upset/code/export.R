# 只封装设备和通用树坐标转换，不隐藏任何统计、布局或配色决策。
options(stringsAsFactors = FALSE)
# 布局预计算也使用支持系统字体的内存设备，避免 pdf(NULL) 的 PostScript 字体警告。
options(device = function(...) ragg::agg_capture(width = 7, height = 7, units = "in", res = 120, background = "white"))

save_gg <- function(plot, stem, width_mm, height_mm, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, device = ragg::agg_png,
                 width = width_mm, height = height_mm, units = "mm", dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, device = grDevices::cairo_pdf,
                 width = width_mm, height = height_mm, units = "mm", bg = "white", family = "Arial")
}

save_draw <- function(expr, stem, width_mm, height_mm, dpi = 300) {
  code <- substitute(expr)
  env <- parent.frame()
  ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm,
                units = "mm", res = dpi, background = "white")
  eval(code, env)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4,
                       height = height_mm / 25.4, family = "Arial", bg = "white")
  eval(code, env)
  grDevices::dev.off()
}

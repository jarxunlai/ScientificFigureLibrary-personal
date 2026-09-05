# =============================================================================
# Nature 2024 Fig. 4a：四物种 GalNAc 通路完整性 UpSet
# organized：线性脚本 + 中文分节
# 数据：论文计数 + 原图列序重建的 exclusive intersection 汇总表
# 不是 UHGG tblastn 复现；作者脚本用 ComplexUpset，本条目用 ggplot2 复刻版式。
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ragg)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
out_dir <- file.path(root, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_font <- "Arial"
species_levels <- c(
  "B. bifidum",
  "C. aerofaciens",
  "F. prausnitzii",
  "F. lactaris"
)
species_cols <- c(
  "B. bifidum" = "#FF0000",
  "C. aerofaciens" = "#2CA25F",
  "F. prausnitzii" = "#F18D00",
  "F. lactaris" = "#4DB6D2"
)
step_levels <- c("Step 0", "Step 1", "Step 2", "Step 3", "Step 4", "Step 5")

# ---- 读入按原图列序整理的 intersection 汇总 --------------------------------
combo <- read_csv(
  file.path(root, "data", "fig4a_intersections.csv"),
  show_col_types = FALSE
)
combo$combo_id <- factor(combo$combo_id, levels = unique(combo$combo_id))
# ggplot2 4.x：因子第一水平在堆积柱顶部。原图自下而上为
# F. lactaris → F. prausnitzii → C. aerofaciens → B. bifidum。
combo$species <- factor(combo$species, levels = species_levels)
step_cols <- c("step0", "step1", "step2", "step3", "step4", "step5")

# ---- set size：每个步骤在多少株中存在 --------------------------------------
set_size <- tibble(
  step = factor(step_levels, levels = step_levels),
  n = vapply(step_cols, function(cn) {
    as.integer(sum(combo[[cn]] * combo$n))
  }, integer(1))
)

# ---- 顶部堆积柱：每个组合按物种计数 ----------------------------------------
bar_df <- combo %>%
  mutate(combo_id = factor(combo_id, levels = unique(combo$combo_id)))

p_bar <- ggplot(bar_df, aes(x = combo_id, y = n, fill = species)) +
  geom_col(width = 0.72, colour = NA) +
  scale_fill_manual(values = species_cols, breaks = species_levels, name = "Species") +
  scale_y_continuous(
    name = "Intersection size",
    breaks = c(0, 1000, 2000, 3000),
    labels = c("0", "1,000", "2,000", "3,000"),
    expand = expansion(mult = c(0, 0.06)),
    limits = c(0, 3800)
  ) +
  theme_classic(base_size = 11, base_family = plot_font) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(8, 8, -6, 2)
  )

p_leg <- cowplot::get_legend(
  p_bar +
    theme(
      legend.position = "left",
      legend.title = element_text(size = 11, family = plot_font),
      legend.text = element_text(face = "italic", size = 10, family = plot_font),
      legend.key.size = unit(0.38, "cm"),
      legend.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(colour = NA)))
)

# ---- 左侧 set-size 水平条 ---------------------------------------------------
p_set <- ggplot(set_size, aes(x = n, y = step)) +
  geom_col(width = 0.62, fill = "#5A5A5A", colour = NA) +
  scale_x_reverse(
    name = "Set size",
    breaks = c(0, 5000, 10000),
    labels = c("0", "5,000", "10,000"),
    limits = c(11000, 0),
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev(step_levels), name = NULL) +
  theme_classic(base_size = 11, base_family = plot_font) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(2, 4, 6, 4)
  )

# ---- 矩阵：灰空心点 + 黑实心点 + 连接线 -------------------------------------
combo_levels <- levels(combo$combo_id)
matrix_df <- combo %>%
  distinct(combo_id, across(all_of(step_cols))) %>%
  pivot_longer(all_of(step_cols), names_to = "step_key", values_to = "present") %>%
  mutate(
    step = factor(
      recode(
        step_key,
        step0 = "Step 0",
        step1 = "Step 1",
        step2 = "Step 2",
        step3 = "Step 3",
        step4 = "Step 4",
        step5 = "Step 5"
      ),
      levels = step_levels
    ),
    combo_id = factor(combo_id, levels = combo_levels),
    present = as.logical(present)
  )

line_df <- matrix_df %>%
  filter(present) %>%
  group_by(combo_id) %>%
  summarise(
    ymin = step[which.min(as.integer(step))],
    ymax = step[which.max(as.integer(step))],
    n_fill = n(),
    .groups = "drop"
  ) %>%
  filter(n_fill > 1)

p_mat <- ggplot(matrix_df, aes(x = combo_id, y = step)) +
  geom_segment(
    data = line_df,
    aes(x = combo_id, xend = combo_id, y = ymin, yend = ymax),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_point(
    data = filter(matrix_df, !present),
    shape = 21,
    size = 3.4,
    fill = "#D9D9D9",
    colour = "#D9D9D9"
  ) +
  geom_point(
    data = filter(matrix_df, present),
    shape = 21,
    size = 3.4,
    fill = "black",
    colour = "black"
  ) +
  scale_y_discrete(limits = rev(step_levels), name = NULL) +
  scale_x_discrete(name = "Step combinations") +
  theme_classic(base_size = 11, base_family = plot_font) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(-6, 8, 6, 0)
  )

# ---- 拼图：左 set-size 对齐矩阵行；图例放在柱图右上空白 --------------------
p <- (wrap_elements(p_leg) + p_bar + plot_layout(widths = c(0.26, 0.74))) /
  (p_set + p_mat + plot_layout(widths = c(0.26, 0.74))) +
  plot_layout(heights = c(1.05, 1))
p <- p +
  plot_annotation(
    title = "a",
    theme = theme(
      plot.title = element_text(
        face = "bold", size = 16, family = plot_font,
        hjust = 0, margin = margin(0, 0, -4, 0)
      ),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

png_path <- file.path(out_dir, "fig4a_galnac_upset.png")
pdf_path <- file.path(out_dir, "fig4a_galnac_upset.pdf")
ragg::agg_png(png_path, width = 10.4, height = 4.6, units = "in", res = 300, background = "white")
print(p)
dev.off()
grDevices::cairo_pdf(pdf_path, width = 10.4, height = 4.6)
print(p)
dev.off()
if (!file.exists(png_path) || file.info(png_path)$size < 1000) {
  stop("PNG 写出失败")
}

file.copy(png_path, file.path(root, "preview.png"), overwrite = TRUE)
cat("PNG", png_path, "\n")
cat("PDF", pdf_path, "\n")
cat("set sizes\n")
print(set_size)
cat("combo totals\n")
print(combo %>% group_by(combo_id) %>% summarise(n = sum(n), .groups = "drop"))

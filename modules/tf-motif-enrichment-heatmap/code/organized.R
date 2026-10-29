# 可移植绘图入口：合成数据示例；不依赖公众号截图或本仓库绝对路径。
# 从模块根目录运行：Rscript --vanilla code/render.R
# 便携模块固定 data/ 相对路径，并写出唯一 canonical 预览 preview.png。
output <- "preview.png"
library(ggplot2)
library(cowplot)
library(grid)
row_meta <- read.csv("data/gene_trends.csv", check.names = FALSE)
motif_meta <- read.csv("data/tf_motifs.csv", check.names = FALSE)
enrichment <- read.csv("data/motif_enrichment.csv", check.names = FALSE)
stopifnot(all(c("row_order", "row_label", "trajectory", "gene_trend") %in% names(row_meta)),
          all(c("motif_order", "motif", "tf_family") %in% names(motif_meta)),
          all(c("row_label", "motif", "adj_p_val", "data_origin") %in% names(enrichment)),
          all(enrichment$data_origin == "independently_synthetic_visual_example"),
          !anyDuplicated(row_meta$row_label), !anyDuplicated(motif_meta$motif),
          all(enrichment$row_label %in% row_meta$row_label),
          all(enrichment$motif %in% motif_meta$motif),
          all(enrichment$adj_p_val >= 0.005 & enrichment$adj_p_val <= 0.045))
row_meta <- row_meta[order(row_meta$row_order), ]
motif_meta <- motif_meta[order(motif_meta$motif_order), ]
row_levels <- rev(row_meta$row_label)
motif_levels <- motif_meta$motif
row_meta$row_label <- factor(row_meta$row_label, levels = row_levels)
motif_meta$motif <- factor(motif_meta$motif, levels = motif_levels)
enrichment$row_label <- factor(enrichment$row_label, levels = row_levels)
enrichment$motif <- factor(enrichment$motif, levels = motif_levels)

# 蓝色方阵和两组独立的类别配色。
family_order <- c("bHLH", "bHLH-ZIP", "bZIP", "C2H2", "ETS", "HD", "HSF",
                  "MADS", "NFAT", "NFI", "NR", "P53L", "RFX", "SMAD", "TBX", "TFAP")
family_colors <- c(bHLH = "#009488", `bHLH-ZIP` = "#FE0000", bZIP = "#13FC00",
                   C2H2 = "#FC2AFD", ETS = "#980000", HD = "#DAFB3A",
                   HSF = "#5FA1CE", MADS = "#D87C5E", NFAT = "#000000",
                   NFI = "#9533FE", NR = "#FEFB04", P53L = "#0037FE",
                   RFX = "#737373", SMAD = "#00FDFE", TBX = "#D0D0D0",
                   TFAP = "#FD9800", unassigned = "#FFFFFF")
trajectory_colors <- c(`L2/3` = "#6D0B14", L4 = "#8B3842", `L5/6` = "#A76A71",
                       CGE = "#1A5700", MGE = "#B2451F", Astro = "#FDC758", ODC = "#265F85")
trend_colors <- c(up = "#F16F48", `trans up` = "#FED231",
                  `trans down` = "#D8D8D5", down = "#00567F")
stopifnot(all(as.character(motif_meta$tf_family) %in% names(family_colors)))

p_matrix_full <- ggplot(enrichment, aes(x = motif, y = row_label, colour = adj_p_val)) +
  geom_point(shape = 15, size = 5, stroke = 0) +
  scale_x_discrete(limits = motif_levels, drop = FALSE, expand = expansion(add = 0.5)) +
  scale_y_discrete(limits = row_levels, drop = FALSE, expand = expansion(add = 0.5)) +
  scale_colour_gradientn(name = "visual value (synthetic)",
                         colours = c("#032F6B", "#77ACD1", "#F2F7FB"),
                         values = scales::rescale(c(0.005, 0.025, 0.045)),
                         limits = c(0.005, 0.045), breaks = c(0.01, 0.02, 0.03, 0.04)) +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0,
                                 barheight = unit(24, "mm"), barwidth = unit(4, "mm"))) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9, base_family = "sans") +
  theme(panel.background = element_rect(fill = "white", colour = NA), panel.border = element_blank(),
        panel.grid.major = element_line(colour = "#A7A7A7", linewidth = 0.20),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length = unit(1.2, "mm"), axis.text = element_blank(),
        legend.position = "right", legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), plot.margin = margin(1.5, 0.5, 0, 0.5, unit = "mm"))
p_value_legend <- suppressWarnings(get_legend(p_matrix_full))
p_matrix <- p_matrix_full + theme(legend.position = "none")

p_left <- ggplot(row_meta, aes(y = row_label)) +
  geom_point(aes(x = "Trajectory", colour = trajectory), shape = 15, size = 4.2, show.legend = FALSE) +
  geom_point(aes(x = "General gene trend", colour = gene_trend), shape = 16, size = 4.8, show.legend = FALSE) +
  scale_x_discrete(limits = c("Trajectory", "General gene trend"), drop = FALSE, expand = expansion(add = 0.5)) +
  scale_y_discrete(limits = row_levels, drop = FALSE, expand = expansion(add = 0.5)) +
  scale_colour_manual(values = c(trajectory_colors, trend_colors)) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_void(base_family = "sans") +
  theme(axis.text.y = element_text(size = 8.8, colour = "black", hjust = 1, margin = margin(r = 2)),
        plot.margin = margin(1.5, 0, 0, 3, unit = "mm"))
annotation_labels <- data.frame(annotation = factor(c("Trajectory", "General gene trend"),
                                                    levels = c("Trajectory", "General gene trend")), y = 0.5)
p_bottom_left <- ggplot(annotation_labels, aes(x = annotation, y = y)) + geom_blank() +
  scale_x_discrete(limits = c("Trajectory", "General gene trend"), drop = FALSE, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 90,
                                   hjust = 1, vjust = 0.5, margin = margin(t = 3)),
        plot.margin = margin(-11, 0, 1, 3, unit = "mm"))

p_family_full <- ggplot(motif_meta, aes(x = motif, y = 0.50, colour = tf_family)) +
  geom_point(shape = 15, size = 4.0) +
  scale_x_discrete(limits = motif_levels, drop = FALSE, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_colour_manual(name = "TF family", values = family_colors, breaks = family_order) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0, ncol = 1,
                               override.aes = list(shape = 15, size = 3.2))) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_void(base_family = "sans") +
  theme(axis.text.x = element_text(size = 8, colour = "black", angle = 90,
                                   hjust = 1, vjust = 0.5, margin = margin(t = 3)),
        legend.position = "right", legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), legend.key.height = unit(3.3, "mm"),
        legend.key.width = unit(3.3, "mm"), plot.margin = margin(-11, 0.5, 1, 0.5, unit = "mm"))
family_legend <- suppressWarnings(get_legend(p_family_full))
p_family <- p_family_full + theme(legend.position = "none")
family_gt <- ggplotGrob(p_family)
family_gt$grobs[grep("axis-b", family_gt$layout$name)] <- lapply(
  family_gt$grobs[grep("axis-b", family_gt$layout$name)],
  function(g) grobTree(g, vp = viewport(y = unit(25, "mm"))))
left_gt <- ggplotGrob(p_bottom_left)
left_gt$grobs[grep("axis-b", left_gt$layout$name)] <- lapply(
  left_gt$grobs[grep("axis-b", left_gt$layout$name)],
  function(g) grobTree(g, vp = viewport(y = unit(25, "mm"))))
aligned_left <- align_plots(p_left, p_bottom_left, align = "v", axis = "lr")
aligned_main <- align_plots(p_matrix, p_family, align = "v", axis = "lr")
left_column <- plot_grid(aligned_left[[1]], left_gt, ncol = 1,
                         rel_heights = c(0.622, 0.378), align = "v", axis = "lr")
matrix_column <- plot_grid(aligned_main[[1]], family_gt, ncol = 1,
                           rel_heights = c(0.622, 0.378), align = "v", axis = "lr")
right_legends <- plot_grid(p_value_legend, family_legend, ncol = 1,
                           rel_heights = c(0.31, 0.69), align = "v", axis = "lr")
p <- plot_grid(left_column, matrix_column, right_legends, nrow = 1,
               rel_widths = c(1.35, 7.65, 1.24), align = "h", axis = "tb")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggsave(output, plot = p, width = 11.39, height = 4.44,
       units = "in", dpi = 300, bg = "white", device = ragg::agg_png)
stopifnot(file.exists(output), file.info(output)$size > 10000)
cat("Rendered synthetic TF motif layout:", normalizePath(output), "\n")

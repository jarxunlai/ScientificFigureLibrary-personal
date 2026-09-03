# =============================================================================
# Nat Commun Fig.5I 背靠背棒棒图（LHSC）
# =============================================================================
# 官方 Source Data Fig5。相对原文：绑定本草稿 xlsx/csv；补 0 行时
# van Galen 组名不再误写成 Pei；x 轴上限改为 20 以免 Pei 超界。
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(ragg)

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else normalizePath(getwd())
  }
)
draft_dir <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
out_dir <- file.path(draft_dir, "output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_font_family <- "Arial"

data <- as.data.frame(read.csv(
  file.path(draft_dir, "data", "fig5.csv"),
  header = TRUE, skip = 1, check.names = FALSE, stringsAsFactors = FALSE
))

genes <- c(
  "CD34", "IL3RA", "SMIM24", "CLEC12A", "IL2RA", "FAM30A", "BEX3", "CD96",
  "CD200", "CDK6", "IL1RAP", "CD33", "SOCS2", "CD9", "KIT", "CD99",
  "CD82", "CPXM1", "CD47", "FCGR2A"
)
data_Leu <- data[data$celltype == "Leukemic Hematopoietic Stem Cell", ]
data_Leu <- data_Leu[data_Leu$marker %in% genes, ]
data_Leu$marker <- factor(data_Leu$marker, levels = rev(genes))

pad_missing <- function(df, group_name) {
  geneno <- setdiff(genes, as.character(df$marker))
  if (!length(geneno)) return(df)
  temp <- data.frame(
    celltype = "Leukemic Hematopoietic Stem Cell",
    group = group_name, marker = geneno, marker_type = "positive",
    gene_impact_score_per_celltype_cell = 0,
    EC_score = 0, specificity = 1
  )
  rbind(df, temp)
}

data_Leu_r <- pad_missing(data_Leu[data_Leu$group == "Pei", ], "Pei")
data_Leu_l <- pad_missing(data_Leu[data_Leu$group == "van Galen", ], "van Galen")
data_Leu_r$marker <- factor(data_Leu_r$marker, levels = rev(genes))
data_Leu_l$marker <- factor(data_Leu_l$marker, levels = rev(genes))

p1 <- ggplot(data_Leu_r, aes(x = gene_impact_score_per_celltype_cell, y = marker)) +
  geom_col(fill = "#a52a2a", width = 0.1) +
  geom_point(aes(size = EC_score), color = "#a52a2a") +
  scale_size_continuous(range = c(0, 8), name = "ECs", breaks = c(2, 4, 6)) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Pei") +
  theme_bw(base_family = plot_font_family) +
  theme(
    plot.title = element_text(face = "bold", color = "#a42929", size = 16),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p2 <- ggplot(data_Leu_l, aes(x = gene_impact_score_per_celltype_cell, y = marker)) +
  geom_col(fill = "#dc968d", width = 0.1) +
  geom_point(aes(size = EC_score), color = "#dc968d") +
  scale_size_continuous(range = c(0, 8), name = "ECs", breaks = c(2, 4, 6)) +
  scale_x_reverse(limits = c(20, 0), breaks = seq(0, 20, 5), expand = expansion(mult = c(0.02, 0))) +
  labs(title = "van Galen") +
  theme_bw(base_family = plot_font_family) +
  theme(
    plot.title = element_text(face = "bold", color = "#a42929", size = 16, hjust = 1),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p <- (p2 | p1) + plot_layout(guides = "collect") &
  theme(
    legend.justification = c("right", "bottom"),
    legend.text = element_text(face = "bold", size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(out_dir, "fig5i_lhsc_lollipop.png"), p, width = 8, height = 5.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(file.path(out_dir, "fig5i_lhsc_lollipop.pdf"), p, width = 8, height = 5.5, device = cairo_pdf, bg = "white")
message("wrote fig5i_lhsc_lollipop.png")

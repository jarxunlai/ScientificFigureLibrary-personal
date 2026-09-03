# =============================================================================
# Cancer Cell 通路 NES 热图
# organized：作者 ComplexHeatmap 代码的可运行版
# =============================================================================
# 与 code/original.R 的可见差异：
#   1. mmc6.xlsx / Colors (ggsci).RData 未随文提供。NES 仍用本草稿合成表
#      （通路名单、星号行号与原文 SelectPathway / InputSig 一致）。
#   2. 注释色用 ggsci pal_npg / pal_cosmic / pal_aaas 近似 ColJournal。
#   3. Heatmap() 参数保持原文：left_annotation、top_annotation、
#      column_split、row_split、rect_gp 白边、cell_fun 星号、
#      height = unit(20,'cm'), width = unit(4,'cm')。
# =============================================================================

library(ComplexHeatmap)
library(circlize)
library(tibble)
library(readr)
library(ggsci)
library(grid)

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

# ---- 第一步：通路名单与 NES 矩阵 -------------------------------------------
SelectPathway <- tibble(
  Pathway = c(
    "KEGG Cell cycle", "KEGG DNA replication", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_FATTY_ACID_METABOLISM", "KEGG Oxidative phosphorylation",
    "HALLMARK_GLYCOLYSIS", "KEGG N-Glycan biosynthesis", "KEGG Lysine degradation",
    "KEGG Apoptosis", "KEGG Autophagy",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS",
    "REACTOME Antigen Presentation", "REACTOME Signaling by Interleukins",
    "REACTOME Cytokine Signaling in Immune system",
    "REACTOME Innate Immune System", "REACTOME Adaptive Immune System",
    "REACTOME Interferon gamma signaling", "HALLMARK_COMPLEMENT",
    "REACTOME Signaling by NOTCH", "REACTOME Signaling by Hedgehog",
    "REACTOME Signaling by FGFR", "KEGG ErbB signaling pathway",
    "KEGG NF-kappa B signaling pathway", "KEGG VEGF signaling pathway",
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING",
    "KEGG Hippo signaling pathway", "REACTOME EPH-Ephrin signaling",
    "REACTOME Signaling by PDGF", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
    "KEGG AMPK signaling pathway", "KEGG HIF-1 signaling pathway",
    "KEGG Tight junction", "REACTOME Cell junction organization",
    "KEGG Focal adhesion",
    "KEGG ECM-receptor interaction",
    "REACTOME Activation of Matrix Metalloproteinases",
    "REACTOME Degradation of the extracellular matrix"
  ),
  Category = c(
    rep("Cell cycle", 3),
    rep("Metabolism", 5),
    rep("Biologic processes", 5),
    rep("Immune response", 7),
    rep("Signalings", 14),
    rep("Cell junction", 3),
    rep("ECM", 3)
  )
)

nes <- as.data.frame(read_csv(file.path(draft_dir, "data", "pathway_nes.csv"), show_col_types = FALSE))
rownames(nes) <- nes$Pathway
Input <- as.matrix(nes[SelectPathway$Pathway, c("A", "B", "C", "D")])
storage.mode(Input) <- "numeric"

# ---- 第二步：星号矩阵（原文 InputSig 行号） --------------------------------
InputSig <- matrix(NA, nrow = nrow(Input), ncol = ncol(Input))
InputSig[c(1:23, 27:37), 1] <- "*"
InputSig[c(1:3, 5, 10:11, 13:26, 28, 30:33, 36, 37), 2] <- "*"
InputSig[c(3, 5:12, 14, 16:19, 25, 27:28, 31, 34:37, 39), 3] <- "*"
InputSig[c(2:5, 8, 11:14, 18, 19, 21, 23:26, 29:32, 37:40), 4] <- "*"

# ---- 第三步：ggsci 近似 ColJournal ------------------------------------------
# Nature NPG: A=F39B7F, B=E64B35, C=3C5488, D=7E6148
ColImmune <- pal_npg("nrc")(10)[c(5, 1, 4, 9)]
names(ColImmune) <- c("A", "B", "C", "D")
cosmic <- pal_cosmic("signature_substitutions")(6)
science8 <- pal_aaas("default")(8)[8]

AnnoCluster <- columnAnnotation(
  Subtype = c("A", "B", "C", "D"),
  col = list(Subtype = ColImmune),
  gp = gpar(col = "white"),
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm")
)
AnnoPathway <- rowAnnotation(
  Category = SelectPathway$Category,
  col = list(
    Category = c(
      "Cell cycle" = cosmic[1],
      "Biologic processes" = cosmic[2],
      "ECM" = cosmic[3],
      "Cell junction" = cosmic[4],
      "Immune response" = cosmic[5],
      "Metabolism" = cosmic[6],
      "Signalings" = science8
    )
  ),
  gp = gpar(col = "white"),
  simple_anno_size = unit(0.5, "cm")
)

# ---- 第四步：原文 Heatmap() -------------------------------------------------
col <- colorRamp2(
  c(-1.5, -1, 0, 1, 1.5),
  c("#00685BFF", "#4CB6ACFF", "white", "#FFDFB2FF", "#EE6C00FF")
)
set.seed(0317)
ht <- Heatmap(
  Input,
  cluster_columns = FALSE,
  name = "NES",
  col = col,
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(9, "cm"),
  top_annotation = AnnoCluster,
  left_annotation = AnnoPathway,
  column_split = c("A", "B", "C", "D"),
  column_title = NULL,
  row_split = SelectPathway$Category,
  row_title = NULL,
  cluster_rows = FALSE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(InputSig[i, j])) {
      grid.text(InputSig[i, j], x, y, gp = gpar(fontsize = 9), just = c("centre", "center"))
    }
  },
  height = unit(40 / 2, "cm"),
  width = unit(4, "cm")
)

png_path <- file.path(out_dir, "pathway_nes_heatmap.png")
pdf_path <- file.path(out_dir, "pathway_nes_heatmap.pdf")
# 原文 pdf 8×12 会让右侧图例压住长通路名；加宽画布并把图例再往右推。
png(png_path, width = 11, height = 12, units = "in", res = 300)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
  gap = unit(5, "mm"),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()
pdf(pdf_path, height = 12, width = 11)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE,
  gap = unit(5, "mm"),
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()
message("wrote ", png_path)

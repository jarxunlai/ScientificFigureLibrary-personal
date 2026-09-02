# =============================================================================
# ggsankeyfier 桑基图：基础 / 渐变 / 多层级 / 代谢通路组合气泡
# organized 版本：中文分节导航，线性脚本，不新增绘图封装函数
# =============================================================================
# 来源：微信公众号《高质量桑基图绘制：布局、配色与组合设计》
#   （R语言数据分析指南，2026-07-27）
#   https://mp.weixin.qq.com/s/PidcbXZv4h56V26IxmkTuA
# 原文未贴完整可运行代码，只给了两段书页截图。本脚本按截图 API
# （ggsankeyfier::pivot_stages_longer / position_sankey / ggplot_build；
#  图 9.60/9.61 用 patchwork 左右拼接）与三张成品图做 visual_inference 复刻。
# 数据为合成示例，不能声称复现书中原始表格或论文结论。
# =============================================================================


library(readr)
library(dplyr)
library(ggplot2)
library(ggsankeyfier)
library(patchwork)
library(RColorBrewer)

# ---- 可移植路径与输出目录 ----------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
data_dir <- file.path(root, "data")
output_arg <- sub("^--output-dir=", "", args[grepl("^--output-dir=", args)])
output_env <- Sys.getenv("SFL_OUTPUT_DIR", unset = "")
out_dir <- if (length(output_arg) && nzchar(output_arg[[1]])) {
  output_arg[[1]]
} else if (nzchar(output_env)) {
  output_env
} else {
  file.path(tempdir(), "ggsankeyfier-layout-color-combo-validation")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
preview_output <- Sys.getenv("SFL_PREVIEW_OUTPUT", unset = "")

sankey_theme <- theme_void(base_family = "") +
  theme(
    plot.margin = margin(6, 18, 10, 6),
    plot.caption = element_text(hjust = 0.5, size = 11, colour = "grey20"),
    legend.position = "none"
  )

# ---- 图 9.60 基础桑基图：脑细胞亚类 -> 预测大类 ---------------------------
# 节点顺序与配色按原文图从像素高度/颜色读取；流量为合成 count。
brain <- read_csv(file.path(data_dir, "brain_subclass_predict.csv"), show_col_types = FALSE)

subclass_levels <- c(
  "IT", "L4 IT", "VIP", "PVALB", "L6 CT", "LAMP5", "SST",
  "Oligodendrocyte", "Astrocyte", "L6b", "L5/6 IT Car3", "L5/6 NP",
  "OPC", "Microglia", "PAX6", "L5 ET", "Endothelial", "Pericyte"
)
predict_levels <- c(
  "Upper-layer intratelencephalic",
  "Deep-layer intratelencephalic",
  "CGE interneuron",
  "MGE interneuron",
  "Miscellaneous",
  "Deep-layer corticothalamic and 6b",
  "Committed oligodendrocyte precursor",
  "Astrocyte",
  "LAMP5-LHX6 and Chandelier",
  "Deep-layer near-projecting",
  "Oligodendrocyte precursor",
  "Microglia",
  "Oligodendrocyte",
  "Vascular"
)

brain_fill <- c(
  "IT" = "#DBB58C",
  "L4 IT" = "#D4B88C",
  "VIP" = "#93C5D9",
  "PVALB" = "#96CAC2",
  "L6 CT" = "#B5C58F",
  "LAMP5" = "#9AC790",
  "SST" = "#96C7D1",
  "Oligodendrocyte" = "#96C8B1",
  "Astrocyte" = "#E6AFA8",
  "L6b" = "#ABC991",
  "L5/6 IT Car3" = "#C6BE8D",
  "L5/6 NP" = "#BEC28E",
  "OPC" = "#97C7A7",
  "Microglia" = "#96C69B",
  "PAX6" = "#C9B48A",
  "L5 ET" = "#A9C39A",
  "Endothelial" = "#8FB9B0",
  "Pericyte" = "#C9B896",
  "Upper-layer intratelencephalic" = "#DBB58C",
  "Deep-layer intratelencephalic" = "#D4B88C",
  "CGE interneuron" = "#93C5D9",
  "MGE interneuron" = "#96CAC2",
  "Miscellaneous" = "#C6BE8D",
  "Deep-layer corticothalamic and 6b" = "#B5C58F",
  "Committed oligodendrocyte precursor" = "#96C8B1",
  "LAMP5-LHX6 and Chandelier" = "#9AC790",
  "Oligodendrocyte precursor" = "#97C7A7",
  "Vascular" = "#8FB9B0"
)

# ggplot y 轴向上；as_is 时因子第一水平在底部。rev() 后与原图自上而下一致。
brain$subclass_label <- factor(brain$subclass_label, levels = rev(subclass_levels))
brain$predict <- factor(brain$predict, levels = rev(predict_levels))
brain <- brain %>% arrange(subclass_label, predict)

brain_long <- pivot_stages_longer(
  brain,
  stages_from = c("subclass_label", "predict"),
  values_from = "count"
)
brain_long$node <- factor(
  as.character(brain_long$node),
  levels = unique(c(rev(subclass_levels), rev(predict_levels)))
)
# 同名节点跨列必须拆开。group 必须是因子：字符 group 会被按字母排序。
# 因子第一水平画在底部，因此水平按原图自上而下的反向排列。
brain_long$group_id <- factor(
  paste(brain_long$stage, brain_long$node, sep = "||"),
  levels = c(
    paste("subclass_label", rev(subclass_levels), sep = "||"),
    paste("predict", rev(predict_levels), sep = "||")
  )
)
brain_long <- brain_long %>% arrange(stage, group_id)
brain_long$fill_col <- unname(brain_fill[as.character(brain_long$node)])
brain_long$label_show <- recode(
  as.character(brain_long$node),
  "Upper-layer intratelencephalic" = "Upper-layer\nintratelencephalic",
  "Deep-layer intratelencephalic" = "Deep-layer\nintratelencephalic",
  "Deep-layer corticothalamic and 6b" = "Deep-layer corticothalamic\nand 6b",
  "Committed oligodendrocyte precursor" = "Committed oligodendrocyte\nprecursor",
  "LAMP5-LHX6 and Chandelier" = "LAMP5-LHX6 and\nChandelier",
  "Deep-layer near-projecting" = "Deep-layer\nnear-projecting",
  "Oligodendrocyte precursor" = "Oligodendrocyte\nprecursor",
  .default = NA_character_
)
brain_long$label_show <- ifelse(
  is.na(brain_long$label_show),
  as.character(brain_long$node),
  brain_long$label_show
)

pos_brain <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.12)
pos_brain_lab_l <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.12, nudge_x = -0.18)
pos_brain_lab_r <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.12, nudge_x = 0.22)

p_brain <- ggplot(
  brain_long,
  aes(x = stage, y = count, group = group_id, connector = connector, edge_id = edge_id)
) +
  geom_sankeyedge(aes(fill = fill_col), position = pos_brain, alpha = 0.55, colour = NA) +
  geom_sankeynode(aes(fill = fill_col), position = pos_brain, colour = NA) +
  geom_text(
    aes(label = ifelse(as.character(stage) == "subclass_label", label_show, NA)),
    stat = "sankeynode",
    position = pos_brain_lab_l,
    hjust = 1,
    size = 2.3,
    family = "",
    lineheight = 0.9,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = ifelse(as.character(stage) == "predict", label_show, NA)),
    stat = "sankeynode",
    position = pos_brain_lab_r,
    hjust = 0,
    size = 2.2,
    family = "",
    lineheight = 0.9,
    na.rm = TRUE
  ) +
  scale_fill_identity() +
  coord_cartesian(clip = "off") +
  labs(caption = "图 9.60　基础桑基图") +
  sankey_theme +
  theme(plot.margin = margin(6, 88, 10, 10))

# ---- 图 9.61 渐变色桑基图：潮间带 -> 生境 ---------------------------------
intertidal <- read_csv(file.path(data_dir, "intertidal_habitat.csv"), show_col_types = FALSE)
hab_levels <- c("Marine", "Freshwater", "Soil", "RefSeq", "Human")
intertidal$habitat <- factor(intertidal$habitat, levels = rev(hab_levels))

it_fill <- c(
  "Intertidal zone" = "#4C9A6A",
  "Marine" = "#C9C6D8",
  "Freshwater" = "#6BAED6",
  "Soil" = "#E8C4A8",
  "RefSeq" = "#B0B0B0",
  "Human" = "#F2D16B"
)

it_long <- pivot_stages_longer(
  intertidal,
  stages_from = c("source", "habitat"),
  values_from = "count"
)
it_long$group_id <- factor(
  paste(it_long$stage, it_long$node, sep = "||"),
  levels = c(
    paste("source", "Intertidal zone", sep = "||"),
    paste("habitat", rev(hab_levels), sep = "||")
  )
)
it_long <- it_long %>% arrange(stage, group_id)
it_long$fill_col <- unname(it_fill[as.character(it_long$node)])

pos_it <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.18)
pos_it_l <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.18, nudge_x = -0.22)
pos_it_r <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.18, nudge_x = 0.22)

p_it <- ggplot(
  it_long,
  aes(x = stage, y = count, group = group_id, connector = connector, edge_id = edge_id)
) +
  geom_sankeyedge(aes(fill = fill_col), position = pos_it, alpha = 0.85, colour = NA) +
  geom_sankeynode(aes(fill = fill_col), position = pos_it, colour = NA) +
  geom_text(
    aes(label = ifelse(as.character(connector) == "from", as.character(node), NA)),
    stat = "sankeynode",
    position = pos_it_l,
    hjust = 1,
    size = 3.2,
    family = "",
    angle = 90,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = ifelse(as.character(connector) == "to", as.character(node), NA)),
    stat = "sankeynode",
    position = pos_it_r,
    hjust = 0,
    size = 3.2,
    family = "",
    na.rm = TRUE
  ) +
  scale_fill_identity() +
  labs(caption = "图 9.61　渐变色桑基图") +
  sankey_theme

p_basic <- p_brain | p_it
ggsave(
  file.path(out_dir, "fig_9_60_9_61.png"),
  p_basic,
  width = 13.2,
  height = 6.8,
  dpi = 300,
  bg = "white"
)

# ---- 图 9.62 多层级桑基图：Global -> Continent -> LandCover -> Habitat -----
# 原文：分层筛选后在基础桑基图上叠加各层级关系。这里用四列宽表一次展开。
globe <- read_csv(file.path(data_dir, "global_landcover_habitat.csv"), show_col_types = FALSE)

continent_levels <- c(
  "North America", "Europe", "Asia", "Australia", "South America",
  "Antarctica", "Africa", "Pacific Ocean", "Atlantic Ocean",
  "Arctic Ocean", "Indian Ocean"
)
landcover_levels <- c(
  "forest", "grassland", "cropland", "aquatic", "desert", "woodland",
  "shrubland", "tundra", "wetland", "urban", "mangrove"
)
habitat_levels <- c(
  "soil", "shoot", "root", "rhizosphere", "deadwood", "air",
  "sediment", "litter", "lichen", "water", "topsoil", "dust"
)

globe$Continent <- factor(globe$Continent, levels = rev(continent_levels))
globe$LandCoverType <- factor(globe$LandCoverType, levels = rev(landcover_levels))
globe$Habitat <- factor(globe$Habitat, levels = rev(habitat_levels))

globe_long <- pivot_stages_longer(
  globe,
  stages_from = c("Global", "Continent", "LandCoverType", "Habitat"),
  values_from = "count"
)
node_levels_globe <- unique(c("Global", rev(continent_levels), rev(landcover_levels), rev(habitat_levels)))
globe_long$node <- factor(as.character(globe_long$node), levels = node_levels_globe)

globe_fill <- c(
  "Global" = "#D07070",
  "North America" = "#7FA86A",
  "Europe" = "#C4A24A",
  "Asia" = "#C47A3A",
  "Australia" = "#9A7A3A",
  "South America" = "#5AA05A",
  "Antarctica" = "#D0C8B8",
  "Africa" = "#C48A4A",
  "Pacific Ocean" = "#6A9EC8",
  "Atlantic Ocean" = "#5B8FBE",
  "Arctic Ocean" = "#8BB8D8",
  "Indian Ocean" = "#4A86B8",
  "forest" = "#5B9A50",
  "grassland" = "#A8C86A",
  "cropland" = "#C4B45A",
  "aquatic" = "#6BB3C8",
  "desert" = "#E0C48A",
  "woodland" = "#7A9A4A",
  "shrubland" = "#C9A05A",
  "tundra" = "#C8D0C0",
  "wetland" = "#6A9A8A",
  "urban" = "#A09090",
  "mangrove" = "#4A8A6A",
  "soil" = "#6A9E78",
  "shoot" = "#8FBF8A",
  "root" = "#7AAB70",
  "rhizosphere" = "#5A9A88",
  "deadwood" = "#A08060",
  "air" = "#B8D4E8",
  "sediment" = "#6A8CA0",
  "litter" = "#C0A878",
  "lichen" = "#A8C0A0",
  "water" = "#6AA0C8",
  "topsoil" = "#8B6B4A",
  "dust" = "#C8B898"
)
globe_long$group_id <- factor(
  paste(globe_long$stage, globe_long$node, sep = "||"),
  levels = c(
    paste("Global", "Global", sep = "||"),
    paste("Continent", rev(continent_levels), sep = "||"),
    paste("LandCoverType", rev(landcover_levels), sep = "||"),
    paste("Habitat", rev(habitat_levels), sep = "||")
  )
)
globe_long <- globe_long %>% arrange(stage, group_id)
globe_long$fill_col <- unname(globe_fill[as.character(globe_long$node)])

pos_g <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.08)
pos_g_lab <- position_sankey(v_space = "auto", order = "as_is", align = "justify", width = 0.08, nudge_x = 0.12)

p_globe <- ggplot(
  globe_long,
  aes(x = stage, y = count, group = group_id, connector = connector, edge_id = edge_id)
) +
  geom_sankeyedge(aes(fill = fill_col), position = pos_g, alpha = 0.45, colour = NA) +
  geom_sankeynode(aes(fill = fill_col), position = pos_g, colour = NA) +
  geom_text(
    aes(label = node),
    stat = "sankeynode",
    position = pos_g_lab,
    hjust = 0,
    size = 2.4,
    family = ""
  ) +
  scale_fill_identity() +
  scale_x_discrete(labels = c("Global", "Continent", "LandCoverType", "Habitat")) +
  labs(caption = "图 9.62　多层级桑基图") +
  sankey_theme +
  theme(
    plot.margin = margin(8, 72, 12, 8),
    axis.text.x = element_text(size = 10, colour = "grey25", margin = margin(t = 4))
  )

ggsave(
  file.path(out_dir, "fig_9_62.png"),
  p_globe,
  width = 9.2,
  height = 7.4,
  dpi = 300,
  bg = "white"
)

# ---- 图 10.24 组合图：代谢物-通路桑基 + 富集气泡 --------------------------
# 书页截图要点：ggplot_build(p1)$data[[4]] 取右侧通路标签 y，
# 再作为气泡图连续 y。截图把 color 写成字符串 "Hit Ratio"；
# 这里按数据列 hit_ratio 映射，并把气泡画进同一坐标系。
links <- read_csv(file.path(data_dir, "metabolite_pathway_links.csv"), show_col_types = FALSE)
enr <- read_csv(file.path(data_dir, "pathway_enrichment.csv"), show_col_types = FALSE)

met_levels <- c(
  "L-Glutamate", "Pyruvate", "2-Oxoglutarate",
  "5-Phospho-alpha-D-ribose 1-diphosphate", "Spermine", "Spermidine",
  "3'-Phosphoadenylyl sulfate", "Adenylyl sulfate", "Deoxyguanosine",
  "GMP", "Xanthosine", "Deoxyinosine", "AMP", "ATP", "GDP",
  "CDP-ethanolamine", "sn-Glycero-3-phosphocholine", "Acetylcholine",
  "Choline", "Choline phosphate", "1-Acyl-sn-glycero-3-phosphocholine",
  "Phosphatidylcholine", "Phosphatidylethanolamine",
  "1-Pyrroline-4-hydroxy-2-carboxylate", "Phosphocreatine", "L-Proline",
  "Citrate", "N-Acetyl-L-aspartate", "Phosphoenolpyruvate", "Isocitrate",
  "Sedoheptulose 7-phosphate", "2-Deoxy-D-ribose 5-phosphate",
  "Galactosylceramide", "N-Acylsphingosine", "Sphingomyelin",
  "Carnosine", "L-Citrulline", "D-Glutamine", "5-L-Glutamyl-taurine",
  "Taurine", "L-Tyrosine"
)
path_levels <- c(
  "Purine metabolism",
  "Glycerophospholipid metabolism",
  "Arginine and proline metabolism",
  "Alanine, aspartate and glutamate metabolism",
  "Citrate cycle (TCA cycle)",
  "Pentose phosphate pathway",
  "Sphingolipid metabolism",
  "beta-Alanine metabolism",
  "Arginine biosynthesis",
  "D-Glutamine and D-glutamate metabolism",
  "Taurine and hypotaurine metabolism",
  "Phenylalanine, tyrosine and tryptophan biosynthesis"
)

met_cols <- c(
  "#55818B", "#67A45E", "#5C9077", "#E4A656", "#6B836B", "#AC864B",
  "#BEA054", "#D0BB5C", "#E5DA66", "#EFEA6F", "#F9E968", "#F4D162",
  "#EBBB5A", "#DE9350", "#517A97", "#517096", "#5C6484", "#6A5871",
  "#794C5E", "#8D424F", "#A63C41", "#C53F39", "#9B5E6E", "#7C548D",
  "#CB7346", "#A9655F", "#689263", "#64A163", "#9F8593", "#AC8298",
  "#B9809E", "#C882A7", "#D382A6", "#BA6F7A", "#8E5C3D", "#598981",
  "#AC6869", "#9F6259", "#7A6A88", "#6A8A7A", "#8A6A5A"
)
names(met_cols) <- met_levels

path_cols <- c(
  "Purine metabolism" = "#C43E39",
  "Glycerophospholipid metabolism" = "#A0C0CE",
  "Arginine and proline metabolism" = "#5B9A50",
  "Alanine, aspartate and glutamate metabolism" = "#DD9190",
  "Citrate cycle (TCA cycle)" = "#ACCD90",
  "Pentose phosphate pathway" = "#F7F2A5",
  "Sphingolipid metabolism" = "#5C4584",
  "beta-Alanine metabolism" = "#B4A4C0",
  "Arginine biosynthesis" = "#E8B47D",
  "D-Glutamine and D-glutamate metabolism" = "#43719D",
  "Taurine and hypotaurine metabolism" = "#D77C47",
  "Phenylalanine, tyrosine and tryptophan biosynthesis" = "#B18772"
)

links$metabolite <- factor(links$metabolite, levels = rev(met_levels))
links$pathway <- factor(links$pathway, levels = rev(path_levels))
links <- links %>% arrange(metabolite, pathway)

combo_long <- pivot_stages_longer(
  links,
  stages_from = c("metabolite", "pathway"),
  values_from = "count"
)
combo_long$node <- factor(
  as.character(combo_long$node),
  levels = c(rev(met_levels), rev(path_levels))
)
combo_long$group_id <- factor(
  paste(combo_long$stage, combo_long$node, sep = "||"),
  levels = c(
    paste("metabolite", rev(met_levels), sep = "||"),
    paste("pathway", rev(path_levels), sep = "||")
  )
)
combo_long <- combo_long %>% arrange(stage, group_id)
combo_long$fill_col <- ifelse(
  combo_long$stage == "metabolite",
  unname(met_cols[as.character(combo_long$node)]),
  unname(path_cols[as.character(combo_long$node)])
)
combo_long$label_show <- recode(
  as.character(combo_long$node),
  "5-Phospho-alpha-D-ribose 1-diphosphate" = "5-Phospho-alpha-D-ribose\n1-diphosphate",
  "1-Acyl-sn-glycero-3-phosphocholine" = "1-Acyl-sn-glycero-\n3-phosphocholine",
  "1-Pyrroline-4-hydroxy-2-carboxylate" = "1-Pyrroline-4-hydroxy-\n2-carboxylate",
  "Glycerophospholipid metabolism" = "Glycerophospholipid\nmetabolism",
  "Arginine and proline metabolism" = "Arginine and proline\nmetabolism",
  "Alanine, aspartate and glutamate metabolism" = "Alanine, aspartate and\nglutamate metabolism",
  "Citrate cycle (TCA cycle)" = "Citrate cycle\n(TCA cycle)",
  "D-Glutamine and D-glutamate metabolism" = "D-Glutamine and\nD-glutamate metabolism",
  "Taurine and hypotaurine metabolism" = "Taurine and hypotaurine\nmetabolism",
  "Phenylalanine, tyrosine and tryptophan biosynthesis" = "Phenylalanine, tyrosine and\ntryptophan biosynthesis",
  .default = NA_character_
)
combo_long$label_show <- ifelse(
  is.na(combo_long$label_show),
  as.character(combo_long$node),
  combo_long$label_show
)

# justify 会把左右两列拉成等高。原图通路终点只约占代谢物柱的一半，
# 气泡框跟通路终点同一段 y。改用 bottom，让右列按节点质量堆在底部。
pos_c <- position_sankey(v_space = "auto", order = "as_is", align = "bottom", width = 0.08)
pos_c_l <- position_sankey(v_space = "auto", order = "as_is", align = "bottom", width = 0.08, nudge_x = -0.16)
pos_c_r <- position_sankey(v_space = "auto", order = "as_is", align = "bottom", width = 0.08, nudge_x = 0.13)

# ggplot2 4.x 会把两个 geom_text 合成一层。第 3 层用独立 geom_rect 占位，
# 第 4 层只放通路标签，对应书页 ggplot_build(p1)$data[[4]]。
p1 <- ggplot(
  combo_long,
  aes(x = stage, y = count, group = group_id, connector = connector, edge_id = edge_id)
) +
  geom_sankeyedge(position = pos_c, fill = "grey80", alpha = 0.55, colour = NA) +
  geom_sankeynode(aes(fill = fill_col), position = pos_c, colour = NA) +
  geom_rect(
    data = data.frame(xmin = 0, xmax = 0, ymin = 0, ymax = 0),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = ifelse(as.character(stage) == "pathway", label_show, NA)),
    stat = "sankeynode",
    position = pos_c_r,
    hjust = 0,
    size = 2.05,
    family = "",
    lineheight = 0.88,
    na.rm = TRUE
  ) +
  scale_fill_identity() +
  scale_x_discrete(expand = expansion(add = c(0.42, 1.35))) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  coord_cartesian(clip = "off") +
  sankey_theme +
  theme(plot.margin = margin(8, 8, 14, 8))

# 书页：lab <- ggplot_build(p1)$data[[4]] %>% select(y, label)
# ggplot2 4.x 可能合并图层，先看第 4 层，没有再按层号往下找。
gb_p1 <- ggplot_build(p1)
layer_order <- seq_along(gb_p1$data)
if (length(gb_p1$data) >= 4) {
  layer_order <- c(4, setdiff(seq_along(gb_p1$data), 4))
}
lab <- NULL
for (layer_i in layer_order) {
  d <- gb_p1$data[[layer_i]]
  if (!("label" %in% names(d))) {
    next
  }
  cand <- d %>%
    filter(!is.na(label), nzchar(as.character(label))) %>%
    mutate(label = gsub("\n", " ", as.character(label))) %>%
    filter(label %in% path_levels) %>%
    select(y, label) %>%
    arrange(desc(y)) %>%
    distinct()
  if (nrow(cand) >= 1) {
    lab <- cand
    break
  }
}
if (is.null(lab)) {
  stop("ggplot_build(p1) 没有取到通路标签 y")
}

# 代谢物标签单独取出后 annotate，避免再占 ggplot_build 层。
met_lab <- ggplot_build(
  ggplot(
    combo_long,
    aes(x = stage, y = count, group = group_id, connector = connector, edge_id = edge_id)
  ) +
    geom_text(
      aes(label = ifelse(as.character(stage) == "metabolite", label_show, NA)),
      stat = "sankeynode",
      position = pos_c_l,
      na.rm = TRUE
    )
)$data[[1]] %>%
  filter(!is.na(label), nzchar(as.character(label)))

# 气泡框贴着通路终点，不拉满整列代谢物高度。
path_nodes <- gb_p1$data[[2]] %>%
  filter(round(x) == 2)
if (nrow(path_nodes) < 1) {
  stop("没有取到通路节点几何")
}
path_span <- max(path_nodes$ymax) - min(path_nodes$ymin)
box_pad <- path_span * 0.02
box_ymin <- min(path_nodes$ymin) - box_pad
box_ymax <- max(path_nodes$ymax) + box_pad

# 气泡画在同一张图里：y 用通路标签坐标，x 把 -log(Pvalue) 映射到通路节点右侧。
logp_min <- 1.0
logp_max <- 5.8
x_box0 <- 2.58
x_box1 <- 3.22
x_scale <- (x_box1 - x_box0 - 0.16) / (logp_max - logp_min)
tick_logp <- 2:5
tick_x <- x_box0 + 0.08 + (tick_logp - logp_min) * x_scale
enr_bub <- enr %>%
  left_join(lab, by = c("pathNames" = "label")) %>%
  mutate(x_bub = x_box0 + 0.08 + (-log(Pvalue) - logp_min) * x_scale)
if (any(is.na(enr_bub$y))) {
  stop("富集表 pathNames 未能匹配 ggplot_build 取出的通路标签")
}

p_combo <- p1 +
  annotate(
    "text",
    x = met_lab$x,
    y = met_lab$y,
    label = met_lab$label,
    hjust = 1,
    size = 2.0,
    family = "",
    lineheight = 0.88
  ) +
  annotate(
    "rect",
    xmin = x_box0,
    xmax = x_box1,
    ymin = box_ymin,
    ymax = box_ymax,
    fill = "white",
    colour = "grey20",
    linewidth = 0.45
  ) +
  geom_point(
    data = enr_bub,
    aes(x = x_bub, y = y, colour = hit_ratio, size = count),
    inherit.aes = FALSE
  ) +
  annotate(
    "segment",
    x = tick_x,
    xend = tick_x,
    y = box_ymin,
    yend = box_ymin + (box_ymax - box_ymin) * 0.018,
    colour = "grey20",
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = tick_x,
    y = box_ymin - (box_ymax - box_ymin) * 0.055,
    label = tick_logp,
    size = 2.3,
    family = "",
    colour = "grey20"
  ) +
  annotate(
    "text",
    x = (x_box0 + x_box1) / 2,
    y = box_ymin - (box_ymax - box_ymin) * 0.12,
    label = "-log(Pvalue)",
    size = 2.6,
    family = "",
    colour = "grey20"
  ) +
  scale_colour_gradientn(
    name = "Hit Ratio",
    colours = colorRampPalette(brewer.pal(9, "YlOrBr")[3:8])(100)
  ) +
  scale_size_continuous(name = "count", range = c(1.8, 6.2)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

ggsave(
  file.path(out_dir, "fig_10_24.png"),
  p_combo,
  width = 14.2,
  height = 8.8,
  dpi = 300,
  bg = "white"
)

# 可选地把生成的组合图写到调用方指定的位置；默认不覆盖输入预览。
if (nzchar(preview_output)) {
  dir.create(dirname(preview_output), recursive = TRUE, showWarnings = FALSE)
  file.copy(
    file.path(out_dir, "fig_10_24.png"),
    preview_output,
    overwrite = TRUE
  )
}

cat("wrote:\n")
cat(" ", file.path(out_dir, "fig_9_60_9_61.png"), "\n")
cat(" ", file.path(out_dir, "fig_9_62.png"), "\n")
cat(" ", file.path(out_dir, "fig_10_24.png"), "\n")

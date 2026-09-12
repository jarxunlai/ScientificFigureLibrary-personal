# 一次性按已验证顺序加载命名空间，然后再附加绘图包。
for (dependency in c("ggtreeExtra", "ggtree", "ade4", "ape", "Biobase", "BiocGenerics", "biomformat", "Biostrings", "cluster", "data.table", "foreach", "igraph", "multtest", "plyr", "reshape2", "scales", "vegan", "phyloseq")) {
  cat("LOAD_NAMESPACE:", dependency, "\n"); flush.console(); loadNamespace(dependency)
}
# 系统发育树与关联注释
# 路径由脚本位置确定；仅使用包内示例，不安装依赖。
# 重要：箱线图对应按 SampleType 合并并按 Order 聚合后的、Abundance < 120 的条件分布。

# 一、依赖和输出位置 -----------------------------------------------------------
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ragg)

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


out_dir <- file.path(root, "output", "figures")
png_path <- file.path(out_dir, "abundance-boxplot-v3.png")
pdf_path <- file.path(out_dir, "abundance-boxplot-v3.pdf")
if (any(file.exists(c(png_path, pdf_path)))) {
  stop("输出已存在；停止以免覆盖。请先审核现有产物并明确选择新的输出名。")
}

# 二、读取包内数据，保留作者原始处理顺序 ---------------------------------------
data("GlobalPatterns")
GP <- GlobalPatterns
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
                              c("Feces", "Skin")
mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order")

# 只展示小于120的观测；不是完整丰度分布，也不是相对丰度百分比。
# 正式解读前须另核对每个类群筛选前后的观测数与有效生物学重复。
melt_simple <- psmelt(mergedGP) %>%
               filter(Abundance < 120) %>%
               # 外圈分类列另命名，避免 geom_fruit 连接树数据时生成 Phylum.x/Phylum.y。
               select(OTU, val=Abundance, abundance_phylum=Phylum)

# 三、只给树附加每个OTU一行的分类信息，避免样本长表重复展开叶端 --------------------
tree <- phy_tree(mergedGP)
tip_phylum <- psmelt(mergedGP) %>%
              distinct(OTU, Phylum) %>%
              rename(label=OTU)
stopifnot(!anyDuplicated(tip_phylum$label),
          setequal(tip_phylum$label, tree$tip.label))
# 只附加叶端分类，不推断内部节点分类，也不改变树分支。
p <- ggtree(tree, layout="fan", open.angle=10) %<+% tip_phylum

# 四、扇形树，叶端与箱体使用同一菌门映射 ---------------------------------------
p <- p +
     geom_tippoint(mapping=aes(color=Phylum),
                   size=1.5,
                   show.legend=FALSE)
p <- rotate_tree(p, -90)

# 五、按树顺序添加外圈箱线图 ---------------------------------------------------
p <- p +
     geom_fruit(
         data=melt_simple,
         geom=geom_boxplot,
         mapping = aes(
                     y=OTU,
                     x=val,
                     group=OTU,
                     fill=abundance_phylum,
                   ),
         size=.2,
         outlier.size=0.5,
         outlier.stroke=0.08,
         outlier.shape=21,
          inherit.aes=FALSE,
         axis.params=list(
                         axis       = "x",
                         text.size  = 1.8,
                         hjust      = 1,
                         vjust      = 0.5,
                         nbreak     = 3,
                     ),
         grid.params=list()
     ) 
     
p <- p +
     scale_fill_discrete(
         name="Phyla",
         guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
     ) +
     theme(
         legend.title=element_text(size=9), 
         legend.text=element_text(size=7) 
     )
p

# 五、显式导出；不自动设为 canonical preview ----------------------------------
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
ggsave(png_path, plot=p, width=9, height=7, units="in", dpi=300,
       bg="white", device=ragg::agg_png)
ggsave(pdf_path, plot=p, width=9, height=7, units="in", bg="white",
       device=grDevices::cairo_pdf)
stopifnot(all(file.info(c(png_path, pdf_path))$size > 0))
packages <- c("ggtreeExtra", "ggtree", "phyloseq", "dplyr", "ggplot2", "ragg")
# 文件生成不等于视觉核验通过；必须检查标签、配色、图例和树-箱线图对齐。

# 已打包的目录自带 preview.png；保留它，并将本次重绘结果写入 output。
preview_path <- file.path(root, "preview.png")
if (!file.exists(preview_path)) {
  stopifnot(file.copy(png_path, preview_path, overwrite=FALSE))
}

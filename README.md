# Open Figure Modules

Open Figure Modules is an openly collaborative collection of portable scientific-figure modules for Scientific Figure Library. Contributions are welcome through pull requests.

The repository keeps cleaned module sources and deterministic ZIP archives together. SFL packages only the derived catalog, previews, thumbnails, and provider notices; complete module ZIP files are fetched only for an exact selected module identity.

## Maintainers

This repository is co-maintained by:

- [jarxunlai](https://github.com/jarxunlai)
- [xuzhougeng](https://github.com/xuzhougeng)

Both maintainers have the same administration, merge, release, and signed-feed publishing authority. Either maintainer may review and merge pull requests, operate maintainer commands, and publish the signed catalog overlay.

GitHub owner-only actions such as deleting or transferring this repository remain with the account that owns it. Day-to-day repository work is shared.

## How Scientific Figure Library uses this repository

- Provider identity: `io.github.jarxunlai.personal-figures`
- Display name: Open Figure Modules
- `main` is the content authority: cleaned sources under `modules/` and generated archives under `archives/`
- `open-figure-feed` is the signed catalog overlay. Installed SFL 0.6.7 and later check `current/source-manifest.json` over HTTPS, verify the Ed25519 signature, and atomically switch the local catalog overlay
- Ordinary module add, update, or withdrawal does not require repackaging the SFL plugin
- Complete module ZIP files are not bundled in the plugin. SFL downloads one exact archive only after the user confirms a module
- Search does not wait for the network. A bundled catalog remains the offline bootstrap until a signed overlay is verified

To use the modules, install Scientific Figure Library and search Open Figure Modules. Do not execute module code from this repository as part of SFL itself.

## Current modules

This repository currently contains **119** modules. The module directories are the public source boundary; each one declares its own metadata, input files, code entry points, preview, and license scope.

当前收录 **119 个绘图模板**。下表按图形类型排列，点击中文名称可查看对应模板的预览、使用说明、示例输入和授权范围。英文名称为可读副标题，而非机器标识。

| Figure family | 模板 / Module | English title |
| --- | --- | --- |
| alluvial | [冲击图](modules/alluvial-plot/) | alluvial plot |
| annotated_heatmap_plus_schematic | [fCpG条形码示意与淋巴肿瘤甲基化热图](modules/nature-evoflux-fig1bc/) | nature evoflux fig1bc |
| back_to_back_lollipop | [NC背靠背棒棒图（LHSC）](modules/nc-backtoback-lollipop/) | Nat Commun back-to-back lollipop plot |
| bar | [柱状抖动散点误差棒](modules/bar-jitter-errorbar/) | bar jitter errorbar |
| bar | [双向柱状图](modules/bidirectional-bar/) | bidirectional bar |
| bar | [复杂百分比柱状图](modules/complex-percent-bar/) | complex percent bar |
| bar | [分面双向柱加热图](modules/facet-bidirectional-bar/) | facet bidirectional bar |
| bar | [渐变柱状图](modules/gradient-bar/) | gradient bar |
| bar | [堆积与百分比柱状图](modules/stacked-percent-bar/) | stacked percent bar |
| benchmark_matrix | [HLCA scIB 整合基准总表](modules/hlca-scib-integration-benchmark-matrix/) | HLCA scIB integration benchmark matrix |
| boxplot | [批量箱线显著性图](modules/batch-boxplot-signif/) | batch boxplot signif |
| boxplot | [箱线图加热图注释](modules/boxplot-heatmap-annot/) | boxplot heatmap annot |
| boxplot | [箱线抖动散点图](modules/boxplot-jitter/) | boxplot jitter |
| boxplot | [带连线箱线图](modules/boxplot-paired-lines/) | boxplot paired lines |
| boxplot | [分半箱线抖动散点图](modules/half-boxplot-jitter/) | half boxplot jitter |
| boxplot | [配对箱线散点连线图](modules/paired-box-scatter-lines/) | paired box scatter lines |
| boxplot | [泛癌箱线显著性图](modules/pancancer-boxplot/) | pancancer boxplot |
| boxplot | [美化箱线图](modules/styled-boxplot/) | styled boxplot |
| bubble | [气泡图加相关性热图](modules/bubble-corr-heatmap/) | bubble corr heatmap |
| bubble | [类别Log2FC气泡图](modules/class-log2fc-bubbles/) | class log2fc bubbles |
| bubble | [复杂气泡图](modules/complex-bubble/) | complex bubble |
| circos | [分类变量环形热图](modules/categorical-circular-heatmap/) | categorical circular heatmap |
| circos | [复杂环形热图](modules/circular-heatmap/) | circular heatmap |
| circos | [环形散点加热图](modules/circular-scatter-heatmap/) | circular scatter heatmap |
| circular_hierarchical_annotation | [海洋微生物 BGC 的 GCC 新颖性与组成概览](modules/nature-2022-fig2a-gcc-overview/) | nature 2022 fig2a gcc overview |
| circular_phylogenomic_tree_with_radial_bars | [海洋微生物 BGC 的系统发育分布](modules/nature-2022-fig2c-phylogenomic-bgc/) | nature 2022 fig2c phylogenomic bgc |
| circular_umap | [环形Circos风格UMAP](modules/umap-style-plot1cell-circlize/) | Circular circos-style UMAP |
| clustered_heatmap_annotation | [差异基因聚类热图加GO和KEGG注释](modules/clustergvis-deg-go-kegg-combo/) | Clustered DEG heatmap with GO and KEGG annotation |
| combo | [堆积柱加折线热图注释](modules/stacked-bar-line-heatmap/) | stacked bar line heatmap |
| deg_heatmap_go_combo | [差异表达热图加GO富集组合图](modules/deg-heatmap-go-panel/) | Combined DEG heatmap and GO enrichment panel |
| density_histogram | [同义替换率密度直方图](modules/ks-density-histogram/) | ks density histogram |
| distribution | [云雨图](modules/raincloud-plot/) | raincloud plot |
| distribution | [雨滴图](modules/raindrop-plot/) | raindrop plot |
| dumbbell | [哑铃图](modules/dumbbell-plot/) | dumbbell plot |
| dumbbell | [多组哑铃图](modules/grouped-dumbbell/) | grouped dumbbell |
| enrichment | [富集分析圈图](modules/enrichment-circle/) | enrichment circle |
| enrichment | [富集气泡图进阶](modules/go-bubble-advanced/) | go bubble advanced |
| enrichment | [KEGG柱状图](modules/kegg-bar/) | kegg bar |
| enrichment_bar | [KEGG变体条形图（蓝黄红柱内标签）](modules/kegg-variant-blue-inside-label/) | KEGG variant bar with blue-yellow-red inside labels |
| enrichment_bar | [KEGG变体条形图（前三白字）](modules/kegg-variant-bw-text-classic/) | KEGG variant bar with white-on-dark labels |
| enrichment_bar | [KEGG变体条形图（柱内通路名）](modules/kegg-variant-inside-label-bar/) | KEGG variant bar with inside-bar pathway labels |
| enrichment_bar_with_genes | [单细胞富集分析条形图（通路+基因）](modules/single-cell-enrichment-bar-pathway-genes/) | Single-cell enrichment bar plot with pathways and genes |
| enrichment_comet_facet | [单细胞GO富集彗星图（分面）](modules/go-enrichment-comet-facet/) | Faceted GO enrichment comet plot |
| enrichment_comet_overlay | [单细胞GO富集彗星图（叠加）](modules/go-enrichment-comet-combined/) | Overlay GO enrichment comet plot |
| faceted_horizontal_bar | [分组水平条形细胞比例图](modules/sc-celltype-grouped-horizontal-bar/) | Faceted horizontal bars of cell-type proportions |
| forest | [分面森林图](modules/facet-forest/) | facet forest |
| forest | [增温效应森林点距图](modules/forest-effect-dotplot/) | forest effect dotplot |
| gantt | [癌症转移生存甘特条](modules/gantt-survival-bars/) | gantt survival bars |
| ggtree_bootstrap_bins | [植物树bootstrap分箱圆点](modules/ggtree-ch13-bootstrap-points/) | ggtree ch13 bootstrap points |
| ggtree_ctldcp_circular | [鸡CTLDcp环形树与禽特异扩张](modules/ggtree-ch13-ctldcp-circular/) | ggtree ch13 ctldcp circular |
| ggtree_fan_antifungal_mutation_tracks | [耳念珠菌系统发育、药敏表型与耐药位点突变](modules/ggtree-ch10-cauris-resistance/) | ggtree ch10 cauris resistance |
| ggtree_genome_locus | [基因组共有基因距离树与位点轨道](modules/ggtree-ch13-genome-locus/) | ggtree ch13 genome locus |
| ggtree_hpv58_distance | [HPV58全基因组树与成对核苷酸距离](modules/ggtree-ch13-hpv58-distance/) | ggtree ch13 hpv58 distance |
| grouped_bar | [分组细胞绝对数 dodge 柱状图](modules/sc-celltype-grouped-dodge-count/) | Grouped dodged bars of absolute cell counts |
| grouped_bar | [样本细胞绝对数 dodge 柱状图](modules/sc-celltype-sample-dodge-count/) | Sample-level dodged bars of absolute cell counts |
| grouped_deg_bar | [Nature风格分组DEG条形图](modules/nature-deg-grouped-barplot/) | Nature-style grouped DEG bar plot |
| heatmap | [Cancer Cell通路NES热图](modules/cancercell-pathway-nes-heatmap/) | Cancer Cell pathway NES heatmap |
| heatmap | [聚类热图加注释](modules/clustered-heatmap-annot/) | clustered heatmap annot |
| heatmap | [复杂热图](modules/complex-heatmap/) | complex heatmap |
| heatmap | [基因融合热图堆积柱](modules/fusion-heatmap-bar/) | fusion heatmap bar |
| heatmap | [热图加渐变连线](modules/heatmap-gradient-links/) | heatmap gradient links |
| heatmap | [热图局部放大](modules/heatmap-local-zoom/) | heatmap local zoom |
| heatmap | [方块热图](modules/tile-heatmap/) | tile heatmap |
| hierarchy | [层状结构可视化](modules/hierarchical-sunburst-tree/) | hierarchical sunburst tree |
| histogram | [分组重叠直方图](modules/overlap-histogram/) | overlap histogram |
| line | [环形分组折线图](modules/circular-grouped-line/) | circular grouped line |
| manhattan | [曼哈顿图](modules/manhattan-plot/) | manhattan plot |
| mantel_corrplot | [Science风格Mantel网络相关热图](modules/science-mantel-corrplot/) | Science-style Mantel network correlation heatmap |
| marker_dotplot | [单细胞marker气泡图加突出框](modules/sc-marker-dotplot-highlight-boxes/) | scRNA marker dotplot with highlight boxes |
| network | [复杂环形互作网络图](modules/circular-interaction-network/) | circular interaction network |
| network | [复杂网络图](modules/complex-network/) | complex network |
| network | [相互作用网络图](modules/interaction-network/) | interaction network |
| open radial multi-track annotation chart | [开放式径向多轨注释图](modules/open-radial-multitrack-annotation/) | open radial multitrack annotation |
| pca | [Nature风格四分组PCA](modules/nature-metabolome-style-pca/) | Nature-style four-group PCA |
| pca | [三维嵌入坐标散点图](modules/pca-3d-scatter/) | pca 3d scatter |
| pca | [PCA图](modules/pca-scatter/) | pca scatter |
| percent_bar | [百分比堆积柱加误差棒](modules/percent-bar-errorbars/) | percent bar errorbars |
| phylogeny | [系统发育树与微生物丰度箱线图](modules/ggtree-ch10-abundance-boxplot/) | ggtree ch10 abundance boxplot |
| pie | [散点饼图](modules/scatterpie/) | scatterpie |
| polar_bar | [南丁格尔玫瑰细胞构成图](modules/sc-celltype-nightingale-rose/) | Nightingale rose chart of cell composition |
| radar | [表达定量雷达图](modules/expression-radar/) | expression radar |
| ranked_nes_scatter | [GSEA打分排序图](modules/hallmark-gsea-nes-ranked-scatter/) | Ranked GSEA NES scatter |
| ridgeline_heatmap | [HOX前后轴山脊图加热图](modules/ncb-fig2d-hox-ridge-heatmap/) | HOX A-P ridgeline with average-expression heatmap |
| sankey | [个性化桑基图](modules/custom-sankey/) | custom sankey |
| sankey | [ggsankeyfier 桑基图布局、配色与组合](modules/ggsankeyfier-layout-color-combo/) | ggsankeyfier Sankey layout, color, and combo |
| sankey | [细胞类型到分组桑基图](modules/sc-celltype-sankey/) | Cell type to group Sankey |
| scatter | [批量散点拟合图](modules/batch-scatter-fit/) | batch scatter fit |
| scatter | [批量散点与折线图](modules/batch-scatter-line/) | batch scatter line |
| scatter | [共定位矩阵散点柱状注释](modules/colocalization-matrix/) | colocalization matrix |
| scatter | [多组学九象限散点图](modules/nine-quadrant-scatter/) | nine quadrant scatter |
| scatter | [散点等高线图](modules/scatter-contour/) | scatter contour |
| scatter | [分组散点直方图注释](modules/scatter-hist-annot/) | scatter hist annot |
| scatter | [复杂形状散点图](modules/shaped-scatter/) | shaped scatter |
| scatter_errorbar | [双向误差棒散点图](modules/scatter-bidirectional-errorbars/) | scatter bidirectional errorbars |
| stacked_area | [跨样本堆叠面积细胞构成图](modules/sc-celltype-stacked-area/) | Stacked area of cell composition across samples |
| stacked_bar | [分组堆叠柱状细胞构成图](modules/sc-celltype-grouped-stacked-bar/) | Grouped 100% stacked bars of cell composition |
| stacked_bar | [样本堆叠比例加总数柱状图](modules/sc-celltype-sample-stacked-proportion/) | Sample stacked proportions with totals |
| stacked_histogram_facets | [海洋微生物 GCF 对参考数据库的新颖性分布](modules/nature-2022-fig2b-gcf-novelty/) | nature 2022 fig2b gcf novelty |
| stacked_proportion_bar | [Nature同款空转生态位堆积柱](modules/nature-spatial-niche-stacked-bar/) | Nature-style spatial niche stacked bar |
| stacked_upset | [四物种GalNAc通路完整性UpSet](modules/nature-fig4a-galnac-upset/) | nature fig4a galnac upset |
| table | [三线表](modules/three-line-table/) | three line table |
| tree | [聚类树](modules/cluster-dendrogram/) | cluster dendrogram |
| umap_density | [暗夜密度热力UMAP](modules/umap-style-density-heatmap/) | Dark magma density UMAP |
| umap_ellipse | [半透明椭圆同色标签UMAP](modules/umap-style-ellipse-labels/) | UMAP with translucent ellipses and matching labels |
| umap_hulls | [单细胞UMAP大群虚线非凸包](modules/umap-unchull-main-type-circles/) | UMAP with dashed non-convex hulls around main cell types |
| umap_in_situ_labels | [簇内直接标注UMAP](modules/umap-style-scrnatoolvis-insitu/) | In-situ labelled UMAP |
| umap_numbered_legend | [编号加侧边图例UMAP](modules/umap-style-scp-numbered-legend/) | Numbered UMAP with count legend |
| umap_square_axes | [方形坐标轴UMAP](modules/umap-style-scp-square-axes/) | Square-axis labelled UMAP |
| umap_stroke | [黑边颗粒风UMAP](modules/umap-style-scpubr-stroke/) | Black-stroke granule UMAP |
| upset | [Upset图](modules/upset-plot/) | upset plot |
| venn | [Venn图](modules/venn-diagram/) | venn diagram |
| violin | [复杂提琴图](modules/complex-violin/) | complex violin |
| violin | [分组提琴显著性图](modules/grouped-violin-signif/) | grouped violin signif |
| volcano | [批量火山图](modules/batch-volcano/) | batch volcano |
| volcano | [气泡火山图](modules/bubble-volcano/) | bubble volcano |
| volcano | [渐变火山图](modules/gradient-volcano/) | gradient volcano |
| volcano | [单细胞分组抖动火山图](modules/grouped-jitter-volcano/) | grouped jitter volcano |
| volcano | [形状火山图](modules/shaped-volcano/) | shaped volcano |
| volcano_go_combo | [单细胞百分比差火山图加GO条形组合图](modules/cell-fig2h-volcano-go-combo/) | Percentage-difference volcano plot with GO bar chart combination |

## Public boundary

Each module declares its own code, content, and documentation licenses. A public repository does not grant redistribution rights to external references. Source screenshots, article or book images, PDFs, real patient or experimental data, unredistributed original code, credentials, and machine-local state are not accepted as module files.

The SFL client only downloads or reads, verifies, extracts, and writes selected modules. It does not run R, Python, notebooks, shell scripts, package installers, or dependency managers.

## Provenance

The modules in the current snapshot were prepared from reviewed Local Published entries. Internal Library revisions, receipts, operation identifiers, locators, and absolute paths are not part of the public repository. Included example inputs and previews demonstrate the plotting layer; they do not claim to reproduce an upstream analysis or scientific conclusion.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md). A pull request is reviewed by an Agent and accepted by a maintainer before it is merged. Either co-maintainer may accept a reviewed pull request. Use a merge commit rather than squash merging: the source commit referenced by each generated archive must remain reachable. Module changes merged into `main` trigger the signed-feed publishing workflow; the updated catalog becomes available only after generation, validation, signing, and publication succeed. Installed SFL clients then verify and adopt the signed overlay independently, without repackaging the plugin. A merged PR alone is not evidence that a client has refreshed its catalog.

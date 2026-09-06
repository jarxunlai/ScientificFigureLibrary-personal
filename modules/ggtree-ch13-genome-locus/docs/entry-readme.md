# 基因组共有基因距离树与位点轨道

## Purpose
问基因内容 Jaccard 距离是否与基因顺序/有无一起变化。

## 使用场景
- 场景一：pangenome / 移动元件比较，树已经分簇，还要看见缺了哪段、多了哪段。
- 场景二：按保守基因对齐轨道，检查插入缺失导致的共线性改变。
- 不适合：MSA 或表达热图。示例数据不能当真实病原体结论。

## Inputs
- `data/example_genes.csv` — gggenes::example_genes

## Commands
`pixi run --environment default Rscript drafts/ggtree-ch13-genome-locus/code/organized.R`

## Limitations
演示数据。确认前不发布。

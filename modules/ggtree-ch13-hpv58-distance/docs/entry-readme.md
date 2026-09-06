# HPV58全基因组树与成对核苷酸距离

## Purpose
问 HPV58 A–D 谱系是否对应全基因组核苷酸距离上的簇。

## 使用场景
- 场景一：全基因组已按命名规则分谱系，要把“树拓扑”和“成对核苷酸距离”放在同一张图里，证明谱系内近、谱系间远。
- 场景二：怀疑某几株分错谱系；右面板距离与颜色不一致时回头查标签/比对。
- 不适合：没有 FASTA、株数极少或只画一棵树。

## Inputs
- `data/HPV58.nwk` — TDbook::tree_HPV58
- `data/HPV58_aln.fas` — TDbook::dna_HPV58_aln

## Methods
groupClade 谱系着色 + geom_facet 成对 Hamming 距离。

## Commands
`pixi run --environment default Rscript drafts/ggtree-ch13-hpv58-distance/code/organized.R`

## Limitations
确认前不发布。

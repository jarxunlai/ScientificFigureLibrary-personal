# HPV58全基因组树与成对核苷酸距离

场景一：一批 HPV / 流感 / SARS-CoV-2 全基因组已按命名规则分成 A/B/C 谱系，审稿人问“谱系是不是只靠树拓扑硬切的”。把树和成对核苷酸距离画在一起：同一谱系的点应聚在左侧、跨谱系拉到右侧，谱系才有序列差异支撑。场景二：怀疑某几株被错分到邻近谱系（重组、污染或 accession 贴错）。右面板里如果某株对“自己谱系”的距离反而比对别的谱系更远，就要回头查标签和比对。不适合：只有 5–6 株、没有比对 FASTA、或只想展示一棵树。距离层需要 N×N 成对距离，株数太多（>200）会糊成一片。

复刻 YuLab treedata-book Figure 13.1 / Chen et al. 2017。树与比对来自 TDbook，未从 GenBank 重下。

This module is an Open Figure Modules submission prepared from a Local Published release. SFL does not execute the code.

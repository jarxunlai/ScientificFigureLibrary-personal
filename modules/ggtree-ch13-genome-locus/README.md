# 基因组共有基因距离树与位点轨道

场景一：比较一组细菌 / 病毒 / 质粒的基因有无（pangenome 或移动元件），Jaccard 距离树已经把基因组分成几簇。想同时看到簇与簇之间到底缺了哪一段、多了哪一段，而不是只看一个距离数字。场景二：怀疑某几株发生了插入/缺失导致共线性改变。把轨道按某个保守基因（本例为 genE）对齐后，就能看出谁丢了上游长基因、谁在中间多了一段。不适合：要展示核苷酸突变或氨基酸保守性（用 MSA）；要展示表达量（用热图）。示例数据是 gggenes::example_genes，不能当成真实病原体结论；换自己的 GFF/基因表时保留 molecule、gene、start、end。

复刻 YuLab treedata-book Figure 13.4。ggtree 4.0.5 geom_motif 会把箭头标成基因组名，已按书中意图改成标 gene 列。

This module is an Open Figure Modules submission prepared from a Local Published release. SFL does not execute the code.

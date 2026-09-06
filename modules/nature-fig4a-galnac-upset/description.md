# 四物种GalNAc通路完整性UpSet

## 需求描述

复刻 Zhernakova et al., Nature 2024 Fig. 4a：UHGG 中四个 ABO 相关肠道菌的 GalNAc 通路步骤完整性 exclusive intersection。问的是这些物种是否真有完整 GalNAc 利用通路，还是只带个别步骤。论文结论：完整通路占 26%（2678/10487）；C. aerofaciens 80% 完整，F. prausnitzii 23% 完整，B. bifidum 几乎没有完整通路。

适用场景（对齐 FigureYa）：场景一——已对 2–6 个物种做完同一条代谢通路的步骤有无，要看完整通路相对残缺组合有多少、缺在哪一步。场景二——宿主基因型与菌种丰度相关，要看是这个种整体都有通路还是只有一部分菌株有；堆积色把物种放进组合柱。不适合韦恩图、基因轨道或生长曲线。

绘图层用 ggplot2 复刻 ComplexUpset 版式。Supplementary Table S16 只有完整通路基因组；长尾组合按原图列序重建，不是 UHGG tblastn 复现。

## 应用场景

未单独记录。此历史模板尚未提供独立应用场景。

## 数据特征

长表：每行 = exclusive intersection × 物种；step0–step5 为 0/1，n 为菌株数。完整通路两行之和必须为 2678。

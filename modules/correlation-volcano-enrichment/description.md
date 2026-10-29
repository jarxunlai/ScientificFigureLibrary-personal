# 相关性火山图与运行富集曲线

Correlation volcanoes and enrichment

## 科学问题

哪些基因或蛋白与连续表型相关，一组通路成员是否集中在相关排序的一端？

Pearson火山图中灰点为背景、青色点高亮WNT成员，FDR阈值虚线与基因标签由输入表映射；下面的运行富集曲线和rug沿同一相关轴对齐。另绘Kendall tau火山图，显著蛋白以暗红突出。

## 应用场景

- 在已有逐基因或逐蛋白相关统计表上高亮关注集合。
- 同时审阅相关方向、显著性与排序富集位置；富集曲线本身不是机制或显著GSEA的证明。

## 输入

合成示例输入；2个运行期输入表。wnt_correlations.csv：1800行×6列；kendall_correlations.csv：1100行×5列

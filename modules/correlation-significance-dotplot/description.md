# 显著性与丰度双编码相关点阵

Correlation dotplots: P value and abundance

## 科学问题

多组特征间的关联方向、统计显著性和丰度如何分别表达？

第一幅按WGD组和通路类别分面，颜色编码Spearman rho、点大小编码-log10P、星号标记显著性；第二幅对细菌门与ARG指标作点阵，颜色编码相关、大小改为相对丰度分档。

## 应用场景

- 比较单细胞转录与基因组改变之间的多层级关联。
- 比较微生物类群与抗性基因指标时，明确区分丰度大小和P值大小两种语义。

## 输入

合成示例输入；2个运行期输入表。wgd_correlations.csv：192行×6列；arg_correlations.csv：150行×5列

# 倒三角相关性热图与聚类树

Inverted triangular correlation heatmap with dendrogram

## 科学问题

哪些变量形成稳定的共变模块，模块之间的关联如何？

**Inverted triangular correlation heatmap with dendrogram**

## 应用场景

代谢物、基因或表型的模块关系概览；聚类树与倒三角热图联合呈现。

输入、绘图步骤与配色均在 `code/organized.R` 中；数据来源和执行边界见 `provenance.json`。

## 输入

合成示例输入；2个运行期输入表。bile_acid_observations.csv：100行×42列；variable_groups.csv：41行×2列

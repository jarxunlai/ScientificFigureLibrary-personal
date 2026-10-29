# 双队列上下三角聚类热图

Dual-cohort triangular clustered correlation heatmap

**Dual-cohort triangular clustered correlation heatmap**

## 科学问题与应用

两个队列的细胞亚群相关结构有哪些一致性与差异？

独立队列或处理组的相关结构对照；双调色板和双聚类树。

输入、绘图步骤与配色均在 `code/organized.R` 中；数据来源和执行边界见 `provenance.json`。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：RColorBrewer, circlize, ggdendro, ggplot2, ragg。

## 输入与执行范围

合成示例输入；2个运行期输入表。GSE165897.csv：120行×29列；In_house.csv：120行×29列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 1 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 04_dual_cohort_clustermap.png

![04_dual_cohort_clustermap.png](preview.png)

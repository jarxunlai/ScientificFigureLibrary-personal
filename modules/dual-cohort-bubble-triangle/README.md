# 上下三角双队列气泡相关性热图

Dual-cohort triangular bubble correlation heatmaps

**Dual-cohort triangular bubble correlation heatmaps**

## 科学问题与应用

M0 与 M1 两类样本的特征相关结构如何在同一坐标框架内比较？

跨分组特征关联；同色版保留分类注释，异色版区分队列与显著性筛选。

输入、绘图步骤与配色均在 `code/organized.R` 中；数据来源和执行边界见 `provenance.json`。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：ComplexHeatmap, RColorBrewer, circlize, ggplot2, ragg。

## 输入与执行范围

合成示例输入；3个运行期输入表。M0_observations.csv：150行×57列；M1_observations.csv：150行×57列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 2 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 09a_same_palette_annotated.png

![09a_same_palette_annotated.png](preview.png)

### 09b_dual_palette_significant.png

![09b_dual_palette_significant.png](figures/09b_dual_palette_significant.png)

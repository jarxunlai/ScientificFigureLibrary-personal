# 差异OTU分类曼哈顿气泡图

Taxonomy-aware differential-OTU Manhattan bubbles

**Taxonomy-aware differential-OTU Manhattan bubbles**

## 科学问题与应用

不同微生物类别在三个比较中的丰度和倍数变化如何分布？

共享科/目排序，气泡编码丰度，空心/实心表示预先给定的显著性。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：dplyr, ggplot2, ggrepel, patchwork, ragg, tidyr。

## 输入与执行范围

合成示例输入；1个运行期输入表。otu.csv：450行×7列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 5 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 17_abundance_legend.png

![17_abundance_legend.png](figures/17_abundance_legend.png)

### 17_B_AZM_vs_B.png

![17_B_AZM_vs_B.png](figures/17_B_AZM_vs_B.png)

### 17_Co_FL_AZM_vs_Co_FL.png

![17_Co_FL_AZM_vs_Co_FL.png](figures/17_Co_FL_AZM_vs_Co_FL.png)

### 17_Co_PS_AZM_vs_Co_PS.png

![17_Co_PS_AZM_vs_Co_PS.png](figures/17_Co_PS_AZM_vs_Co_PS.png)

### 17_otu_manhattan.png

![17_otu_manhattan.png](preview.png)

# 调控基因集放射网络与外圈数量

Radial regulon network with outer counts

*Radial regulon network with outer counts*

## 科学问题与应用

不同功能类别下的基因集规模及层级归属如何同时展示？

转录调控基因集分类；功能模块层级与规模展示，不把层级关系解读为实测调控边。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：ragg。

## 输入与执行范围

合成示例输入；1个运行期输入表。gene_sets.csv：68行×3列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 2 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 05_regulon_network_bars.png

![05_regulon_network_bars.png](preview.png)

### 05_regulon_network_bubbles.png

![05_regulon_network_bubbles.png](figures/05_regulon_network_bubbles.png)

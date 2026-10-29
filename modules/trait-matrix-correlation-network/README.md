# 分组性状散点矩阵与相关网络

Grouped trait matrix and network

下三角按森林类型颜色和叶习性形状展示五个叶片性状的散点与拟合，右上方网络用红蓝区分负正相关、实虚线区分显著与不显著，并以线宽表达相关强度。

## 科学问题与应用

叶片性状的分布模式与生理指标相关网络是否支持同一观察？

- 联合查看生态性状的观测分布和性状—生理指标关联。
- 跨生态系统展示关系时，同时保留分组符号以识别组间混杂。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：dplyr, ggplot2, patchwork, ragg, tidyr。

## 输入与执行范围

合成示例输入；1个运行期输入表。leaf_traits.csv：95行×10列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 1 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 07_grouped_scatter_matrix_network.png

![07_grouped_scatter_matrix_network.png](preview.png)

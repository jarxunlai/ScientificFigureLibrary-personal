# 倒三角相关热图与聚类树

Inverted triangle and dendrogram

将相关矩阵的上三角以菱形单元旋转为倒三角，顶部保留complete-linkage聚类树，两条斜边分别标注相同变量，并用类别色点提示所属分组。

## 科学问题与应用

变量是否形成相近的共变模块，模块边界如何与分类对应？

- 展示共表达、代谢物或生理指标的模块化相关结构。
- 用较紧凑的倒三角版式保留树与变量名称，供后续人工拼版。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：dplyr, ggdendro, ggplot2, patchwork, ragg, scales, tidyr。

## 输入与执行范围

合成示例输入；2个运行期输入表。abundance.csv：90行×27列；groups.csv：26行×2列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 1 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 05_inverted_triangle_dendrogram.png

![05_inverted_triangle_dendrogram.png](preview.png)

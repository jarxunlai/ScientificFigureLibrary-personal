# 对数散点边缘直方图与残差箱线图

Scatter, marginals and residuals

主图以双对数坐标显示基因组大小与CDS数，按SAG、MAG、WGS着色并叠加log尺度线性拟合；上方和右侧是同源数据的边缘直方图，图内残差箱线图显示三组的拟合偏离。

## 科学问题与应用

基因组长度与编码序列数的关系及其残差是否因数据类型不同而变化？

- 比较不同组装来源的基因组尺度关系。
- 对两变量拟合后检查分组残差，而不把相关关系解释为因果。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：dplyr, ggplot2, patchwork, ragg, scales, tidyr。

## 输入与执行范围

合成示例输入；1个运行期输入表。genomes.csv：1180行×3列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 1 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 06_scatter_marginal_hist_residual.png

![06_scatter_marginal_hist_residual.png](preview.png)

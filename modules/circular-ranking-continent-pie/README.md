# 国家排名外圈与洲别贡献饼图

Country ranking bars with a continent contribution pie

*Country ranking bars with a continent contribution pie*

## 科学问题与应用

国家贡献排名与洲别总体占比如何保持数量一致？

全球贡献结构；国家-洲别两层空间汇总和排名。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：ragg。

## 输入与执行范围

合成示例输入；1个运行期输入表。country_contributions.csv：90行×5列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 4 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 10_country_continent_source_blue.png

![10_country_continent_source_blue.png](figures/10_country_continent_source_blue.png)

### 10_country_continent_source_green.png

![10_country_continent_source_green.png](figures/10_country_continent_source_green.png)

### 10_country_continent_source_orange.png

![10_country_continent_source_orange.png](figures/10_country_continent_source_orange.png)

### 10_country_continent_tutorial_purple.png

![10_country_continent_tutorial_purple.png](preview.png)

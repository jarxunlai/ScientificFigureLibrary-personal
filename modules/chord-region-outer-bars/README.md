# 作物区域和弦与管理措施外圈条形图

Crop-region chord diagram with management bars

*Crop-region chord diagram with management bars*

## 科学问题与应用

作物、区域与驱动因子的关系如何与管理措施指标并置？

农业管理措施概览；两套类别体系下的关系与指标联合展示。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：ragg。

## 输入与执行范围

合成示例输入；3个运行期输入表。inner_links.csv：162行×3列；crop_regions.csv：18行×3列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 2 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 06_crop_region_chord.png

![06_crop_region_chord.png](figures/06_crop_region_chord.png)

### 06_crop_region_chord_management.png

![06_crop_region_chord_management.png](preview.png)

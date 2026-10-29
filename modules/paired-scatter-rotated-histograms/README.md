# 配对散点与对角线旋转直方图

Paired scatter with diagonal rotated histograms

**Paired scatter with diagonal rotated histograms**

## 科学问题与应用

处理前后配对读数的差异在条件和动物之间如何变化？

两只动物三条件的六个面板；显示配对检验和差值中位数。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：cowplot, dplyr, ggplot2, ggrepel, patchwork, ragg, tidyr。

## 输入与执行范围

合成示例输入；1个运行期输入表。paired_neural.csv：330行×6列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 7 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 18_Monkey_D_Combined.png

![18_Monkey_D_Combined.png](figures/18_Monkey_D_Combined.png)

### 18_Monkey_D_Vestibular.png

![18_Monkey_D_Vestibular.png](figures/18_Monkey_D_Vestibular.png)

### 18_Monkey_D_Visual.png

![18_Monkey_D_Visual.png](figures/18_Monkey_D_Visual.png)

### 18_Monkey_F_Combined.png

![18_Monkey_F_Combined.png](figures/18_Monkey_F_Combined.png)

### 18_Monkey_F_Vestibular.png

![18_Monkey_F_Vestibular.png](figures/18_Monkey_F_Vestibular.png)

### 18_Monkey_F_Visual.png

![18_Monkey_F_Visual.png](figures/18_Monkey_F_Visual.png)

### 18_neural_paired_six_panel.png

![18_neural_paired_six_panel.png](preview.png)

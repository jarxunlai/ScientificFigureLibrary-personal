# 分组散点与线性对数拟合六面板

Grouped scatter with linear and logarithmic fits

**Grouped scatter with linear and logarithmic fits**

## 科学问题与应用

不同生态系统中初始存量与损失、干扰时间与相对损失有什么关联？

三种生态系统、四类干扰；分别绘制线性和log(x)回归。

## 运行

```text
Rscript --vanilla code/organized.R <output-directory>
```

在已经安装依赖的 R 环境中运行；代码不会安装包、访问网络或运行原研究上游流程。输入从本模块定位，不依赖当前工作目录。输出目录通过第一个参数明确传入；不传入时使用 R 会话临时目录。图形维持模板布局与配色；Arial 字体可用时保持该字体。

依赖：dplyr, ggplot2, ggrepel, patchwork, ragg, tidyr。

## 输入与执行范围

合成示例输入；1个运行期输入表。blue_carbon.csv：300行×7列。每个输入的列名和行数见 data_schema.yml。公开版已实际运行，绘图执行范围为合成数据；不是原论文数值或实验结果。

## 独立图像

本模板对应 7 幅独立 PNG，同时生成配对 PDF；它们不合并成一张替代预览。PDF 为运行输出，不随本公共源码分发。

### 12_blue_carbon_six_panel.png

![12_blue_carbon_six_panel.png](preview.png)

### 12_Mangrove_absolute.png

![12_Mangrove_absolute.png](figures/12_Mangrove_absolute.png)

### 12_Mangrove_relative.png

![12_Mangrove_relative.png](figures/12_Mangrove_relative.png)

### 12_Saltmarsh_absolute.png

![12_Saltmarsh_absolute.png](figures/12_Saltmarsh_absolute.png)

### 12_Saltmarsh_relative.png

![12_Saltmarsh_relative.png](figures/12_Saltmarsh_relative.png)

### 12_Seagrass_absolute.png

![12_Seagrass_absolute.png](figures/12_Seagrass_absolute.png)

### 12_Seagrass_relative.png

![12_Seagrass_relative.png](figures/12_Seagrass_relative.png)

# 双向误差棒散点图

## Purpose

同时显示 x、y 方向标记分数区间，比较是否落在 y=x 附近。

## Visual structure

散点加横向和纵向误差棒，虚线斜率为 1，形状区分基因型。

## Data requirements

xmin/xmax/ymin/ymax 由 runif 模拟。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 模拟数据；白色点在白底上几乎看不见，这是原配色。

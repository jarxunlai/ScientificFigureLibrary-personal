# 同义替换率密度直方图

## Purpose

叠加多对比较的 Ks 分布，比较同义替换率峰值。

## Visual structure

直方图加密度曲线，多组半透明填充叠在同一横坐标。

## Data requirements

脚本内 rnorm 模拟 8 列 Ks；无外部表。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 原脚本第一层直方图引用了未定义的 tmp_data，整理版改为 after_stat(density)。
- 数据为模拟，不是真实比对得到的 Ks。

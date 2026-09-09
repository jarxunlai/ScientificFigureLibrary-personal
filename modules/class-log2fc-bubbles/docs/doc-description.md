# 类别Log2FC气泡图

## Purpose

按脂质类别展示 log2FC 与显著性，点大小表示 -log10(P)。

## Visual structure

横向气泡图，垂直虚线在 0，点填充按类别着色。

## Data requirements

宽表模拟后成长表 Class / Log2FC / Pvalue。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 数据和 p 值由 runif 模拟，每次运行不同。

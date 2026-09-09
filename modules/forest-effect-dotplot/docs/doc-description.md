# 增温效应森林点距图

## Purpose

按类群展示增温效应大小、区间和显著性。

## Visual structure

翻转坐标的点距图；背景色带分组；颜色区分正负效应与显著性。

## Data requirements

rownames.csv 提供 28 个类群名；效应值和 p 标记由脚本随机生成。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- min/med/max 与显著性符号每次运行都会变，预览只代表一次随机实现。
- 不是原文统计模型复现。

# 癌症转移生存甘特条

## Purpose

用甘特条比较四种癌的原发灶与转移灶 40 个月总生存曲线下面积。

## Visual structure

横向甘特条按癌种分行，色块表示转移类型，点为均值，右侧标注样本数。

## Data requirements

每行一种癌×转移类型；low/high 为区间，metastasis 为类型，number 为样本数。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 示例数据来自课程包 test_data.csv，不是原文生存分析复现。
- 原脚本用 ggplot2 默认 linewidth 前的 size 参数画磁贴边框。

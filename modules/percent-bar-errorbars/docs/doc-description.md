# 百分比堆积柱加误差棒

## Purpose

展示多样本祖先成分比例，并加上转换后的误差棒。

## Visual structure

翻转的 100% 堆积柱，灰/黄/蓝三组分，黑色误差棒。

## Data requirements

30 个样本 × 3 组分由 sample() 生成。

## Recommended use

本地查看该图类型的 ggplot2 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 均值和标准差每次运行都会变。
- 误差棒是课程示例算法，不是真实 ADMIXTURE 标准误。

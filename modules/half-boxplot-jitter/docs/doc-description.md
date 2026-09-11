# 分半箱线抖动散点图

## Purpose

复现 KS 生信绘图课程包「041分半箱线图+抖动散点图」的绘图层。

## Visual structure

见 `preview.png` 与课程原图；未做出版级视觉审查。

## Recommended use

本地查看该图类型的 ggplot2/R 实现；需要同类布局时可改数据后重跑 `code/organized.R`。

## Limitations

- 示例数据来自课程包或脚本内模拟，不是原文分析复现。
- 本地 Gallery 收录授权仅覆盖课程图代码与布局参考；科学结论与出版级视觉审查未做。
- gghalves 的 transformation= 在 ggplot2 4 下失败，整理版用箱线+jitter 近似。

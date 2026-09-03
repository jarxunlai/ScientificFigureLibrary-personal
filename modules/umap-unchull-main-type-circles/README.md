# 单细胞UMAP大群虚线非凸包

**UMAP with dashed non-convex hulls around main cell types**

A UMAP module with colored main-type hull outlines, subcluster labels, and a compact coordinate cue.

## 用途

Show broad cell-type territories around subclusters in a synthetic embedding.

## 输入

Synthetic cell-level coordinates with subcluster and main-cell-type columns.

## 运行

需要 R 4.5 或兼容版本，以及：ggplot2, dplyr, ggrepel, ggunchull, tidydr, grid。

```text
Rscript code/example.R
```

脚本只读取本模块的相对路径示例数据，并把生成结果写到模块输出目录或临时目录；不会访问外部数据、Local Published 内部状态或本机绝对路径。

## 公开边界

本模块只包含清洗后的代码、示例数据、生成预览和文档。来源截图、参考图片、PDF、日志、内部 receipt、未确认再分发的原始代码以及本机状态不在模块内。示例数据用于展示绘图层，不代表真实患者、实验样本、论文数据或科研结论。

SFL materialize 只下载、校验、解压和写入锁文件；不会运行 R、Python、shell、notebook、安装器或依赖管理器。

# 单细胞GO富集彗星图（叠加）

**Overlay GO enrichment comet plot**

A portable overlay comet plot for pathway enrichment rows across multiple cell clusters.

## 用途

Compare enrichment strength and contributing-gene counts for several cell clusters on a shared axis.

## 输入

A fully synthetic long table with Cluster, Description, Count, pvalue, and -log10(p) fields.

## 运行

需要 R 4.5 或兼容版本，以及：dplyr, readr, ggplot2, ggforce, ragg。

```text
Rscript code/example.R
```

脚本只读取本模块的相对路径示例数据，并把生成结果写到模块输出目录或临时目录；不会访问外部数据、Local Published 内部状态或本机绝对路径。

## 公开边界

本模块只包含清洗后的代码、示例数据、生成预览和文档。来源截图、参考图片、PDF、日志、内部 receipt、未确认再分发的原始代码以及本机状态不在模块内。示例数据用于展示绘图层，不代表真实患者、实验样本、论文数据或科研结论。

SFL materialize 只下载、校验、解压和写入锁文件；不会运行 R、Python、shell、notebook、安装器或依赖管理器。

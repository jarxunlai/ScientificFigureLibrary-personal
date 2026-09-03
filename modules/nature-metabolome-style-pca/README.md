# Nature风格四分组PCA

**Nature-style four-group PCA**

A four-group PCA presentation with confidence ellipses, outlined points, and direct sample labels.

## 用途

Show group separation and sample labels in a clean PCA panel.

## 输入

Synthetic sample-level coordinates with sample, group, PC1, and PC2 columns.

## 运行

需要 R 4.5 或兼容版本，以及：ggplot2, ggrepel, readr, ragg。

```text
Rscript code/example.R
```

脚本只读取本模块的相对路径示例数据，并把生成结果写到模块输出目录或临时目录；不会访问外部数据、Local Published 内部状态或本机绝对路径。

## 公开边界

本模块只包含清洗后的代码、示例数据、生成预览和文档。来源截图、参考图片、PDF、日志、内部 receipt、未确认再分发的原始代码以及本机状态不在模块内。示例数据用于展示绘图层，不代表真实患者、实验样本、论文数据或科研结论。

SFL materialize 只下载、校验、解压和写入锁文件；不会运行 R、Python、shell、notebook、安装器或依赖管理器。

# Nature Fig.2b — GCF 新颖性分布

本目录是独立 Gallery 条目：`nature-2022-fig2b-gcf-novelty`。未来 Local Published 必须单独发布，不与 Fig.2a/c 合并。

## Purpose

比较 GCF 到 RefSeq/BiG-FAM 和 MIBiG 的平均最小 cosine distance，并按 BGC class 和 phylum 解读新颖性来源。

## Inputs

- `data/gcf-distance-summary.csv`：6,907 个独立 GCF 距离。
- `data/histogram-observations.csv`：GCF-category assignment 长表。
- 原始来源：官方 Supplementary Table 2，见 `provenance.yml`。

## Methods

作者兼容的非加权 GCF×category 堆叠直方图，30 bins，黑线为 dbar=0.2。独立计数复算为 RefSeq-novel 3,861、MIBiG-novel 6,688；y 轴写为 GCF-category count，避免把混合类别的堆叠高度误当独立 GCF 数。

## Software and data sources

实际 R 会话见 `validation/render-session.txt`。已观察：R 4.5.3；ggplot2 4.0.3；dplyr 1.2.1；tidyr 1.3.2；patchwork 1.3.2；ragg 1.5.2。Wisp/model exact version unavailable。

## Commands and scripts

```powershell
pixi run --locked --environment default Rscript --vanilla gallery/nature-2022-fig2b-gcf-novelty/code/organized.R
```

从项目根目录运行。脚本不依赖 panel a 或 c，图例独立。

## Outputs

- `preview.png`
- `output/figures/fig2b-gcf-novelty.png`
- `output/figures/fig2b-gcf-novelty.pdf`
- `data/histogram-bins.csv`
- `validation/paper-number-checks.csv`

## Limitations

公开表/代码类别数和正文有 3000/3012、1816/1815 两处差异；未改数据凑数。未运行 Local Published 导入/发布，未 Git commit。

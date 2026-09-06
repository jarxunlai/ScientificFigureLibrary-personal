# 鸡CTLDcp环形树与禽特异扩张

## Purpose
问鸡 CTLDcp 是否在 Group V 形成禽特异基因扩张。

## 使用场景
- 场景一：跨物种基因家族树，看某物种拷贝是一对一直系同源还是亚家族扩张。
- 场景二：同时标 Group II/V 和 Subgroup；tip 多时用环形省版面。
- 不适合：一两个物种或时间树。tip 前缀必须能区分物种。

## Inputs
- `data/CTLDcp.nwk` — TDbook::tree_treenwk_30.4.19

## Commands
`pixi run --environment default Rscript drafts/ggtree-ch13-ctldcp-circular/code/organized.R`

## Limitations
确认前不发布。

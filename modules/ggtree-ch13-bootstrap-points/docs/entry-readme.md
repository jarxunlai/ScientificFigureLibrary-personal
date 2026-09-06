# 植物树bootstrap分箱圆点

## Purpose
标出哪些内部分叉有高/中/低 bootstrap 支持。

## 使用场景
- 场景一：带 bootstrap 的 Newick 要在论文里标明哪些科/属节点可靠；三档圆点比全写数字更干净。
- 场景二：近缘种是否单系——连接节点若是黑点（BP<70），不能写成定论。
- 不适合：没有 node.label；后验概率更适合连续色条。

## Inputs
- `data/RMI_tree.nwk` — TDbook::text_RMI_tree

## Commands
`pixi run --environment default Rscript drafts/ggtree-ch13-bootstrap-points/code/organized.R`

## Limitations
确认前不发布。

# fCpG条形码示意与淋巴肿瘤甲基化热图

## 需求描述

复刻 Gabbutt et al., Nature 2025 Fig. 1b–c：用示意棒棒糖解释 bulk 甲基化分布如何编码群体历史，再用 978×2,204 热图展示淋巴 fCpG 的斑驳条形码。作者仓库无 Fig. 1b/c 脚本；本条目是 Source Data + 期刊图的视觉复刻，不是 EVOFLUx 推断复现。

## 应用场景

- 场景一：只有 bulk 450k / EPIC 甲基化矩阵，要向读者解释为什么中间 β 值的单峰不等于“半甲基化位点”，而是多克隆平均。把 Fig. 1b 上排（星形树 + 杂乱灰点 + 单峰）和下排（W 形）并排放，比单独贴直方图更不容易被读成批次效应。
- 场景二：已经筛出一批候选 fCpG，需要证明它们不像 clock CpG 那样按病种成块，而是跨病种斑驳。Fig. 1c 用官方列序 + 顶部 Sample group / tumour fraction / discovery / cancer / platform 五条注释，适合作为“位点是条形码、不是细胞类型标志”的主图。
- 场景三：论文要同时交代发现集（Training=True）和验证/健康对照（Training=False）。顶部青/橙色条能直接标出哪些列参与了 fCpG 发现，避免把健康对照误当成发现集。
- 不适合：要展示差异甲基化（DMP/DMH）或按病种聚类的表观亚型；那种图应让位点按疾病成块，而不是斑驳。也不适合单细胞甲基化 UMAP，或只有几十个样本的小热图。

## 数据特征

官方 Source Data Fig. 1c 宽表（注释行 + fCpG 行）；Fig. 1b 三列模拟 bulk 分数；棒棒糖长表。fCpG 矩阵不进 Git，运行时从 inbox xlsx 读取。

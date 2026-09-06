# fCpG条形码示意与淋巴肿瘤甲基化热图

复刻 [Gabbutt, Duran-Ferrer et al., *Nature* 2025](https://doi.org/10.1038/s41586-025-09374-4) **Fig. 1b–c**。

根据 PubMed（PMID 40931062，PMC12443617）正文与期刊图注，这两块回答的问题是：

> bulk 甲基化阵列上的 fluctuating CpG（fCpG）能不能当“自然条形码”，把多克隆、近期克隆扩张、扩张后继续波动三种群体历史区分开？在 2,204 个淋巴样本上，这些位点是按病种成块，还是每个肿瘤各自斑驳？

答案是：**能区分，而且是斑驳条形码，不是病种标志。** 多克隆群体把 0 / 0.5 / 1 三种等位基因状态平均成约 0.5 的单峰；近期克隆扩张把创始细胞的状态写成 W 形（峰在 0、0.5、1）；扩张后再波动会把 W 形抹宽。跨样本热图因此呈 speckled pattern，健康淋巴对照则是均匀中间甲基化。这支持后续用 EVOFLUx 从单次 bulk 剖面推断进化历史，但 Fig. 1b–c 本身还不给出生长速率或预后。

## 生物学问题（读图时抓住什么）

1. **Fig. 1b**：同一套 fCpG，为什么 bulk 直方图会长成单峰或 W 形？原因是细胞谱系结构，不是探针失败。
2. **Fig. 1c**：978 个淋巴 fCpG 在 2,204 个样本上是否按 B-ALL / CLL / MM 成块？不成块，说明位点编码的是谱系噪声，可作跨病种条形码。
3. **放大框**：肿瘤块红蓝相间；健康块接近白/浅粉，对应中间 β。这是发现标准的阳性对照，不是装饰。

## 图的构成

**Panel b（示意，三行）**

| 行 | 树 | 单细胞棒棒糖 | bulk 直方图 |
| --- | --- | --- | --- |
| Polyclonal fCpGs | 星形，远 MRCA | 每列一个细胞，灰白黑不规则 | 单峰，中心 ~0.5 |
| Clonal fCpGs | 长干 + 短末端 | 同行全同，Mean methylation 为 0 / 0.5 / 1 | W 形 |
| Evolving fCpGs | 层次树 | 大部分锁定，少数位点已漂 | 仍可见 0 / 0.5 / 1，峰更宽 |

白 = 双等位基因未甲基化，灰 = 杂合，黑 = 双等位基因甲基化。

**Panel c（数据热图）**

1. 主体：978 行 fCpG × 2,204 列样本，蓝低红高，中心 0.5 为白。
2. 顶部注释（自上而下）：Sample group（17 类）、Tumour fraction、fCpG discovery（True/False）、Sample type（Cancer/Normal）、Platform（450k/EPIC）。
3. 列序来自官方 Source Data，已是层次聚类结果。
4. 右侧放大：Tumour samples 取中右斑驳区（列 1280–1460）；Healthy samples 取右侧连续健康块（列 2135–2172）。

## 适用场景

写法对齐本仓库已发布条目（如 ggtree bootstrap 分箱）和 FigureYa 模块说明：先写“什么实验/数据已经在手”，再写“这张图让读者先看懂哪一句生物学判断”，最后写“不要拿它去回答什么”。

- **场景一：方法文或综述要解释“中间甲基化 ≠ 半甲基化位点”。** 手里只有 bulk β 矩阵时，单独贴直方图常被读成批次或纯度问题。把 Fig. 1b 三行树 + 棒棒糖 + 直方图并排放，读者能把单峰/W 形直接映射到谱系结构。
- **场景二：已经筛出候选 fCpG / 随机波动位点，要证明它们不是细胞类型或病种标志。** 用官方列序的大热图：若按病种成块，位点更像 clock 或 lineage marker；若跨 B-ALL、CLL、MM 斑驳，才支持“条形码”叙事。顶部 Sample group 色条是读图锚点。
- **场景三：发现集和验证集混在一张图里。** fCpG discovery 青/橙色条标出 Training=True 的 1,471 列；健康对照和部分复发/缓解样本为 False。适合在方法图里预先交代哪些列进入了位点筛选。
- **场景四：需要同时交代阵列平台和肿瘤纯度。** 450k（紫）与 EPIC（绿）混排、纯度从白到深蓝，用来回答“斑驳是不是平台或低纯度假象”。期刊图上平台并不按病种成块。

## 不适用

- 差异甲基化热图、表观亚型热图、或按基因通路聚类的甲基化图：那些图必须让位点按生物学分组，而不是斑驳。
- 只有几十个样本、或只有健康对照：W 形和 speckled pattern 都看不出来。
- 单细胞甲基化 UMAP / haplotype 图：那是 Extended Data Fig. 2 的问题，不是 Fig. 1c。
- 声称复现 EVOFLUx 推断、生长速率或 CLL 预后：那些是 Fig. 2–5。

## 数据与代码来源

- 官方 Source Data Fig. 1：`41586_2025_9374_MOESM5_ESM.xlsx`（sheet `Figure 1b` / `Figure 1c`）。
- 样本元数据：`CalumGabbutt/evoflux` 的 `data/BloodMethMetadata.csv`（2,204 行，与 Fig. 1c 列完全重合）。
- 作者代码：无 Fig. 1b/c 脚本。Python 包只做推断；Duran-Ferrer 仓库只覆盖 Fig. 1g、4、5。
- 完整 β 矩阵另存 Zenodo `10.5281/zenodo.15479736`；本条目不下载、不重跑发现流程。

## 状态

- `code/original.R` 记录“作者无此图脚本”。
- `code/organized.R` 绑定官方 Source Data。
- 仅 `private_reference`。确认前不晋升 gallery，也不把期刊原图当 canonical preview。

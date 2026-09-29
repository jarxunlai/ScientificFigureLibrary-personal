# 大鼠心血管细胞 marker–GO 组合点图（Supplementary Figure 2A 视觉复刻）

本模块展示 26 个细胞簇 × 130 个 marker 的点阵和与之对齐的 GO 术语三联 NES 条形。灵感和细胞类型/基因/GO 标签来源：Arduini et al., *Cell Reports* 44, 115091 (2025), DOI [10.1016/j.celrep.2024.115091](https://doi.org/10.1016/j.celrep.2024.115091)，Supplementary Figure 2A。独立生成的预览已去掉多面板字母 `a`。

**仅视觉示范，不是作者分析或真实实验值。** 点的表达比例、表达色值及 GO NES 全部由固定种子 `20250128` 生成（`origin=synthetic_visual_mimic`）；不可据此论证细胞表达或通路显著性。论文原图、PDF、Table S2、原始矩阵均不随模块分发。论文全文标称 CC BY-NC-ND 4.0；本模块含来源图中转录的标签和视觉编排。公开再分发权限由项目所有者在本次任务中声明具有完整授权；wisp-science 未核验授权文件，此声明不更改原论文自身许可证。独立 R 代码以 MIT 许可发布；拟态 CSV、标签示例和新生成的预览拟按 CC BY 4.0 发布，前提是所有者对上述再分发与许可拥有相应权限。

## 输入与绘图方法

- `data/panel_a_labels.csv`：26 簇的显示名、取样颜色、每簇 5 个基因和 3 个 GO 术语；文字转录来自论文图。
- `data/dotplot_synthetic.csv`：26×130 点阵；点大小为合成百分比，颜色为按基因 0–1 缩放的合成表达。
- `data/go_nes_synthetic.csv`：每簇 3 条合成 NES，按显示顺序与每 5 个基因行对齐；不代表真实 fgsea 结果。

`code/plot.R` 按行展开 5 个 marker；基因作为纵轴，黑色簇名纵排在点阵上方，中列独立彩色簇名，每组三条条形长度映射 NES（横轴 0–3），GO 名称置于条框之外。两框均有轮廓。画幅 249×350 mm，PNG 300 dpi；绘图调用 R base graphics 与 `ragg::agg_png`。

## 运行和输出

从**模块根目录**执行 `Rscript --vanilla code/plot.R`，需要 R 和 ragg。输出 `output/preview.png` 与 `output/preview.pdf`。项目 Windows 环境中的实际测试：`pixi run --locked --environment default Rscript --vanilla code/plot.R`（从模块根运行）；R 4.5.3、ragg 1.5.2。其他 OS 和版本没有运行验证。

## 来源、权利与局限

- 源文献 DOI 见上；原论文图仅为本地参考，**不是**本模块的预览输入文件。
- 真实 GEO `GSE280111` 原始数据 >1 GB，本模块未运行作者的差异分析、表达量统计或 GO GSEA；作者未提供原始绘图代码。
- `LICENSE` 仅适用于本模块自写 R 代码；`DATA_LICENSE` 说明其余分发资产的提议许可与来源限制。再分发声明基于用户陈述，而非 wisp-science 的独立版权审查。

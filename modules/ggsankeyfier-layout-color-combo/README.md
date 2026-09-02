# ggsankeyfier 桑基图布局、配色与组合

这是一个可移植的个人 R 绘图模块，展示四类版式：两列基础桑基图、带间渐变桑基图、四列多层级桑基图，以及代谢物–通路桑基与富集气泡组合图。

## 公开边界

- `code/organized.R` 是根据公开可见版式和 API 线索整理的线性示例脚本，关系标记为 `visual_inference`，不是原作者完整代码。
- `data/` 中的五张表是为版式演示构造的合成数据，不代表原始研究数据、真实样本或论文结论。
- `preview.png` 和 `thumbnail.jpg` 是本模块清洗内容生成的预览图。
- 本模块的 MIT / CC BY 4.0 许可只覆盖本仓库中明确列出的清洗代码、合成数据、文档和生成预览；外部文章、书页、参考图片、PDF、第三方字节和未公开原始代码不在分发边界内。

来源线索：<https://mp.weixin.qq.com/s/PidcbXZv4h56V26IxmkTuA>。来源页面没有提供完整可运行脚本和原始表格，因此本模块只记录可审查的版式参考，不声称复现来源数据。

## 运行

需要 R 4.5 或兼容版本，以及 `ggplot2`、`ggsankeyfier`、`dplyr`、`readr`、`tidyr`、`patchwork` 和 `RColorBrewer`。

```text
Rscript code/organized.R
```

脚本从自身位置解析 `data/`，默认把验证 PNG 写入 R 的临时目录，不会在模块树中生成未声明文件。也可以通过 `SFL_OUTPUT_DIR` 或 `--output-dir=...` 指定输出目录，并通过 `SFL_PREVIEW_OUTPUT` 指定生成的组合预览输出位置。脚本不会执行安装器，也不会依赖某台机器的绝对路径。

## 输入和输出

- `data/brain_subclass_predict.csv`：脑亚类到预测大类的合成流向。
- `data/intertidal_habitat.csv`：潮间带到生境的合成流向。
- `data/global_landcover_habitat.csv`：Global–Continent–LandCoverType–Habitat 合成路径。
- `data/metabolite_pathway_links.csv`：代谢物到通路的合成关联。
- `data/pathway_enrichment.csv`：通路富集气泡的合成数值。
- `preview.png`：主预览；`thumbnail.jpg`：搜索缩略图。

## 状态和限制

发布者验证只覆盖合成数据上的脚本可执行性和版式检查；它不等于对外部科研结果、原始数据或科学结论的再验证。SFL materialize 阶段只读取、校验和解压固定 ZIP，不运行本模块代码。

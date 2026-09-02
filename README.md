# ScientificFigureLibrary-personal

这是 Scientific Figure Library 的个人内容仓库。源码模块和确定性 ZIP 共用同一仓库；SFL 插件只携带派生 Catalog、预览和缩略图，不携带完整 ZIP。

## 当前模块

- `modules/ggsankeyfier-layout-color-combo/`：经清洗的 ggsankeyfier 桑基图示例，包含合成 CSV、整理版 R 代码、主预览和缩略图。
- `archives/ggsankeyfier-layout-color-combo.zip`：由 source commit 的模块文件树确定性生成的归档。
- `catalog/`：归档清单和逐模块准入记录。

模块级许可证、发布者审核状态和公开边界以各模块 `module.yml` 与 README 为准；不要把仓库根目录的 MIT 代码许可证自动套用到模块内容。

## 内容边界

本仓库不包含书页、文章截图、PDF、参考图片、真实患者/实验数据、未授权原始代码或本机环境状态。SFL materialize 只校验和解压固定归档，不执行模块代码。

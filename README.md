# Open Figure Modules

Open Figure Modules is an openly collaborative collection of portable scientific-figure modules for Scientific Figure Library. Contributions are welcome through pull requests.

The repository keeps cleaned module sources and deterministic ZIP archives together. SFL packages only the derived catalog, previews, thumbnails, and provider notices; complete module ZIP files are fetched only for an exact selected module identity.

## Current modules

This repository currently contains **36** modules. The module directories are the public source boundary; each one declares its own metadata, input files, code entry points, preview, and license scope.

- `cancercell-pathway-nes-heatmap`
- `cell-fig2h-volcano-go-combo`
- `clustergvis-deg-go-kegg-combo`
- `deg-heatmap-go-panel`
- `ggsankeyfier-layout-color-combo`
- `go-enrichment-comet-combined`
- `go-enrichment-comet-facet`
- `hallmark-gsea-nes-ranked-scatter`
- `hlca-scib-integration-benchmark-matrix`
- `kegg-variant-blue-inside-label`
- `kegg-variant-bw-text-classic`
- `kegg-variant-inside-label-bar`
- `nature-deg-grouped-barplot`
- `nature-metabolome-style-pca`
- `nature-spatial-niche-stacked-bar`
- `nc-backtoback-lollipop`
- `ncb-fig2d-hox-ridge-heatmap`
- `sc-celltype-grouped-dodge-count`
- `sc-celltype-grouped-horizontal-bar`
- `sc-celltype-grouped-stacked-bar`
- `sc-celltype-nightingale-rose`
- `sc-celltype-sample-dodge-count`
- `sc-celltype-sample-stacked-proportion`
- `sc-celltype-sankey`
- `sc-celltype-stacked-area`
- `sc-marker-dotplot-highlight-boxes`
- `science-mantel-corrplot`
- `single-cell-enrichment-bar-pathway-genes`
- `umap-style-density-heatmap`
- `umap-style-ellipse-labels`
- `umap-style-plot1cell-circlize`
- `umap-style-scp-numbered-legend`
- `umap-style-scp-square-axes`
- `umap-style-scpubr-stroke`
- `umap-style-scrnatoolvis-insitu`
- `umap-unchull-main-type-circles`

## Public boundary

Each module declares its own code, content, and documentation licenses. A public repository does not grant redistribution rights to external references. Source screenshots, article or book images, PDFs, real patient or experimental data, unredistributed original code, credentials, and machine-local state are not accepted as module files.

The SFL client only downloads or reads, verifies, extracts, and writes selected modules. It does not run R, Python, notebooks, shell scripts, package installers, or dependency managers.

## Provenance

The modules in the current snapshot were prepared from reviewed Local Published entries. Internal Library revisions, receipts, operation identifiers, locators, and absolute paths are not part of the public repository. Included example inputs and previews demonstrate the plotting layer; they do not claim to reproduce an upstream analysis or scientific conclusion.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md). A pull request is reviewed by an Agent and accepted by the maintainer before it is merged. Merge alone does not automatically place content in an SFL package; deterministic archives and the bundled snapshot are generated and verified separately.

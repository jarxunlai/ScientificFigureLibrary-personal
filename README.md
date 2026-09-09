# Open Figure Modules

Open Figure Modules is an openly collaborative collection of portable scientific-figure modules for Scientific Figure Library. Contributions are welcome through pull requests.

The repository keeps cleaned module sources and deterministic ZIP archives together. SFL packages only the derived catalog, previews, thumbnails, and provider notices; complete module ZIP files are fetched only for an exact selected module identity.

## Maintainers

This repository is co-maintained by:

- [jarxunlai](https://github.com/jarxunlai)
- [xuzhougeng](https://github.com/xuzhougeng)

Both maintainers have the same administration, merge, release, and signed-feed publishing authority. Either maintainer may review and merge pull requests, operate maintainer commands, and publish the signed catalog overlay.

GitHub owner-only actions such as deleting or transferring this repository remain with the account that owns it. Day-to-day repository work is shared.

## How Scientific Figure Library uses this repository

- Provider identity: `io.github.jarxunlai.personal-figures`
- Display name: Open Figure Modules
- `main` is the content authority: cleaned sources under `modules/` and generated archives under `archives/`
- `open-figure-feed` is the signed catalog overlay. Installed SFL 0.6.7 and later check `current/source-manifest.json` over HTTPS, verify the Ed25519 signature, and atomically switch the local catalog overlay
- Ordinary module add, update, or withdrawal does not require repackaging the SFL plugin
- Complete module ZIP files are not bundled in the plugin. SFL downloads one exact archive only after the user confirms a module
- Search does not wait for the network. A bundled catalog remains the offline bootstrap until a signed overlay is verified

To use the modules, install Scientific Figure Library and search Open Figure Modules. Do not execute module code from this repository as part of SFL itself.

## Current modules

This repository currently contains **46** modules. The module directories are the public source boundary; each one declares its own metadata, input files, code entry points, preview, and license scope.

- `cancercell-pathway-nes-heatmap`
- `cell-fig2h-volcano-go-combo`
- `clustergvis-deg-go-kegg-combo`
- `deg-heatmap-go-panel`
- `ggsankeyfier-layout-color-combo`
- `ggtree-ch13-bootstrap-points`
- `ggtree-ch13-ctldcp-circular`
- `ggtree-ch13-genome-locus`
- `ggtree-ch13-hpv58-distance`
- `go-enrichment-comet-combined`
- `go-enrichment-comet-facet`
- `hallmark-gsea-nes-ranked-scatter`
- `hlca-scib-integration-benchmark-matrix`
- `kegg-variant-blue-inside-label`
- `kegg-variant-bw-text-classic`
- `kegg-variant-inside-label-bar`
- `nature-2022-fig2a-gcc-overview`
- `nature-2022-fig2b-gcf-novelty`
- `nature-2022-fig2c-phylogenomic-bgc`
- `nature-deg-grouped-barplot`
- `nature-evoflux-fig1bc`
- `nature-fig4a-galnac-upset`
- `nature-metabolome-style-pca`
- `nature-spatial-niche-stacked-bar`
- `nc-backtoback-lollipop`
- `ncb-fig2d-hox-ridge-heatmap`
- `open-radial-multitrack-annotation`
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

See [`CONTRIBUTING.md`](CONTRIBUTING.md). A pull request is reviewed by an Agent and accepted by a maintainer before it is merged. Either co-maintainer may accept a reviewed pull request. Merge alone does not automatically place content in an SFL package; deterministic archives and the signed catalog overlay are generated and verified separately.
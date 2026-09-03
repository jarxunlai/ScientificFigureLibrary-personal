# Open Figure Modules

Open Figure Modules is an openly collaborative collection of portable scientific-figure modules for Scientific Figure Library. Contributions are welcome through pull requests.

The repository keeps cleaned module sources and deterministic ZIP archives together. SFL packages only the derived catalog, previews, thumbnails, and provider notices; complete module ZIP files are fetched only for an exact selected module identity.

## Current modules

- `ggsankeyfier-layout-color-combo`
- `sc-celltype-grouped-horizontal-bar`
- `sc-celltype-sankey`
- `sc-celltype-grouped-stacked-bar`
- `sc-celltype-stacked-area`
- `sc-celltype-grouped-dodge-count`
- `sc-celltype-sample-dodge-count`
- `sc-celltype-sample-stacked-proportion`
- `sc-celltype-nightingale-rose`

## Public boundary

Each module declares its own code, content, and documentation licenses. A public repository does not grant redistribution rights to external references. Source screenshots, article or book images, PDFs, real patient or experimental data, unredistributed original code, credentials, and machine-local state are not accepted as module files.

The SFL client only downloads or reads, verifies, extracts, and writes selected modules. It does not run R, Python, notebooks, shell scripts, package installers, or dependency managers.

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md). A pull request is reviewed by an Agent and accepted by the maintainer before it is merged. Merge alone does not automatically place content in an SFL package; deterministic archives and the bundled snapshot are generated and verified separately.

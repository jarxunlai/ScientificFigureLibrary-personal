# Maintainer notes

Generated archive and SFL snapshot commands are run from the SFL core checkout:

```text
npm run modules:validate -- --check --repository <this-repository>
npm run modules:archive -- --write --repository <this-repository>
npm run modules:catalog -- --write --repository <this-repository>
npm run modules:source-pack -- --write --repository <this-repository>
```

The commands do not create commits, push, run R, install dependencies, or modify the Gallery. Review generated diffs before each Git operation.

## CI

`open-figure-feed-validate` does **not** run `modules:archive --check`.

That check rebuilds every ZIP from HEAD, rewrites `catalog/archive-manifest.json`, and compares ZIP bytes. Archives submitted by ScientificFigureLibrary on Windows cannot be byte-reproduced on Ubuntu UTC (ZIP DOS timestamps). Pull-request validation instead:

1. `validate-personal-modules --check` for licenses, portable paths, and excluded private files
2. `publish-open-figure-feed --mode validate`, which unzips committed archives and checks them against each module's recorded `sourceCommit`

Use `modules:archive --write` locally only when you intend to regenerate archives.

# Maintainer notes

Generated archive and SFL snapshot commands are run from the SFL core checkout:

```text
npm run modules:validate -- --check --repository <this-repository>
npm run modules:archive -- --write --repository <this-repository>
npm run modules:catalog -- --write --repository <this-repository>
npm run modules:source-pack -- --write --repository <this-repository>
```

The commands do not create commits, push, run R, install dependencies, or modify the Gallery. Review generated diffs before each Git operation.

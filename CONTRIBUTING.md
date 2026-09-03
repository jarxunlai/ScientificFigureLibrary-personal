# Contributing to Open Figure Modules

Contributions are welcome through GitHub pull requests.

## Module source changes

- Add or update only a complete `modules/<moduleId>/` source tree.
- Use a stable, readable, lowercase module ID.
- Declare the exact file inventory, `requiredFiles`, code entry point, inputs, preview identity, thumbnail identity, licenses, publisher review, and execution scope in `module.yml`.
- Keep plotting examples portable: resolve the module root from the script location, use relative inputs, and write outputs only to an explicit output directory or a temporary directory.
- Do not run installers from module code.

## Files that must not be submitted

- source or reference screenshots, article images, paper figures, book pages, or PDFs without explicit redistribution permission;
- patient data, real experimental data, private matrices, credentials, tokens, keys, or machine-local paths;
- `.pixi/`, R session state, caches, environment logs, validation scratch files, receipts, or Local Library state;
- original third-party code without a clear compatible redistribution license.

## Licenses and provenance

Every module must identify licenses separately for code, content, and documentation. By submitting a pull request, the contributor confirms that they have the right to contribute the submitted bytes under the licenses declared in the module manifest. External material may be cited textually, but its bytes must not be included unless redistribution permission is established.

Contributor attribution is preserved through the pull request and Git history. Personal email addresses or internal identity fields are not added to the runtime catalog.

## Review and generated files

The maintainer or an assigned Agent validates paths, file inventory, manifests, preview identities, public boundaries, and deterministic archive generation. The maintainer makes the final merge decision.

Do not hand-edit files under `archives/` or the aggregate archive manifest. Those files are generated from an exact source commit after source review. Pull requests use two retained logical commits: a source commit followed by a generated archive commit. Do not squash these commits when merging.

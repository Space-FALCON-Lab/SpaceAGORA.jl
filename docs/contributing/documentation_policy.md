# Documentation Policy

How SpaceAGORA's documentation is organized, what keeps it in step with the
code, and which checks a change must pass. To run SpaceAGORA, see the
[user guide](../src/index.md); for the rest of the contributor material, see
the [contributor guide](index.md).

Only user workflows and the supported API and extension references belong in
`docs/src/`, the source of the published site. Contributor procedures and
repository contracts stay outside it. The documentation falls into three kinds,
each kept current in its own way.

## Generated API docs

The API reference pages are generated from Julia docstrings and from the list of
public names in:

```text
docs/public_api_symbols.jl
```

These pages are the reference for the exported package surface. A pull request
that changes the root exports updates `docs/public_api_symbols.jl` in the same
pull request; `test/gates/ci_public_api_surface_gate.jl` fails when the list and
the exports differ. A listed name without a docstring fails the docs build (see
[documentation enforcement](#documentation-enforcement)).

## Contract docs

The architecture and quality contracts in `docs/architecture/` and
`docs/quality/` are hand-written and normative. They state the repository's
rules, ownership boundaries and quality gates, and the generated API pages do
not replace them.

## Operational docs

Operational docs describe workflows: the CLI, assets and licensing boundaries,
distributed and HPC launches, extension recipes and templates. Each should be
tied either to a test gate or to a machine-readable source of truth, as these
already are:

- Assets: `data/assets_manifest.toml`, which the asset catalog and the CLI read.
- The generated API pages: `docs/public_api_symbols.jl`.
- The extension hook reference: generated from `docs/extension_hooks.jl`, which
  `test/gates/ci_extension_hook_reference_gate.jl` checks against the package.
- The distributed/HPC, extensibility and API policy pages, with the templates,
  starters and HPC scripts they describe:
  `test/gates/ci_hpc_extensibility_docs_gate.jl` checks that those files exist
  and that the pages keep their key terms.
- The worked examples in `docs/src/user/gnc_feedback.md` and
  `docs/src/user/gnc_planning.md`, and the two extension walkthroughs
  (`custom_effector.md` and `custom_actuator.md`): every region of a script in
  `examples/gnc_algorithm_env/` or `examples/walkthrough_env/` between a
  `# snippet: <name>` line and `# end snippet` must appear on the script's page
  as one Julia code block, character for character. CI runs the scripts, and
  `test/gates/ci_hpc_extensibility_docs_gate.jl` compares them with the pages,
  so change the page and the script together.

Terminal commands in user pages should work on every platform, or come with a
variant for each shell where the syntax differs.

## PR review checklist

Every pull request goes through the docs drift items in
`.github/pull_request_template.md`:

- if root exports changed, update `docs/public_api_symbols.jl`;
- if assets changed, update `data/assets_manifest.toml` and `docs/src/assets.md`;
- if new extension points were added, update `docs/src/extensibility.md`;
- if a hook was added, renamed or deprecated, update `docs/extension_hooks.jl`;
- if new CLI commands or flags were added, update `docs/src/cli.md`.

## Release gate

Before merge to `main`, the documentation release gate is:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
julia --project=. test/gates/ci_public_api_surface_gate.jl
julia --project=. test/gates/ci_hpc_extensibility_docs_gate.jl
```

Treat a failure of any of them as blocking. CI runs the four in the
`documentation-release-gate` job of `.github/workflows/julia-ci.yml`, where the
build gets `--no-deploy` so that it publishes nothing. The contract
suite, `test/contracts/pr_runtests.jl`, runs the two gates again, together with
`test/gates/ci_extension_hook_reference_gate.jl` and
`test/gates/ci_docs_contract_gate.jl`. The last checks that the pull request
template and this page keep the items above, that `docs/make.jl` keeps the
site's navigation, and that contributor pages stay out of the published site.
The Documentation workflow, `.github/workflows/docs.yml`, builds the site again
and publishes it from `main`.

## Documentation enforcement

The docs build fails when a name listed in `docs/public_api_symbols.jl` has no
docstring. The only exceptions are the names in `KNOWN_PLACEHOLDER_DOCS` in the
same file, and the build also fails when one of those has since been documented
or is no longer public. The build runs Documenter with `checkdocs = :exports`
and `warnonly = false`, so any Documenter warning, such as an unresolved
`@ref`, fails it too.

The build writes its report of undocumented names to:

```text
docs/build/undocumented_public_exports.txt
```

The report covers only the names in `docs/public_api_symbols.jl`, not internal
modules. The Documentation workflow uploads it as an artifact of each run.

## Build the docs locally

From the repository root, use Julia 1.12 to prepare the docs environment and
build the documentation. The setup may download packages:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
```

The docs environment, `docs/Project.toml`, holds Documenter and takes
SpaceAGORA from this checkout; the root project does not include Documenter.
On every run, `docs/setup.jl` copies the root manifest to the ignored
`docs/Manifest.toml` and resolves from it, so the build uses the package
versions the tests use. Run it again after the root manifest changes.
`docs/make.jl` stops with an error when it is started from another environment,
before the setup has run, or when SpaceAGORA loads from a different checkout.

The build checks external links, and GitHub limits anonymous requests, so a
local link check can fail with HTTP 429. Give the build a GitHub token: in bash
or zsh, `GITHUB_TOKEN=$(gh auth token) julia --project=docs docs/make.jl`; in
PowerShell, set `$env:GITHUB_TOKEN = gh auth token` first. Documenter sends the
token only with its requests to GitHub.

After a successful build, open `docs/build/index.html`. Pages under
`docs/contributing/`, `docs/architecture/` and `docs/quality/` are not part of
that site; when you change one, read its Markdown, its commands and its
repository links directly.

Generated documentation stays local, in the ignored `docs/build/`, `docs/site/`
and `docs/src/generated/` directories.

## Figures

Figures live in `docs/src/assets/figures/` and are committed; the
documentation build does not regenerate them. The scripts in
`scripts/docs_figures/` reproduce them, and that directory's README lists, for
each figure, the page that uses it and what its script needs. When a change
alters what a figure shows, rerun its script and commit the result. Diagrams
are SVG, plots PNG and viewer screenshots JPEG, which keeps the files small.
Every image needs alt text that says what it shows.

## Regenerate common local artifacts

One helper recreates the common ignored reports and the documentation:

```text
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```

Each mode runs its study or build, so none is needed to inspect simulation
results that already exist. `runtime-analysis` prepares the examples
environment first and writes to `output/performance/` unless given another
directory, and it and `telemetry` also accept `full` in place of `quick`. The
`docs` mode runs the two documentation commands above.

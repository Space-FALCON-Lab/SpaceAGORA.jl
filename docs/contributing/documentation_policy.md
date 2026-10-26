# Documentation Policy

Use this page when you are changing the docs site, the generated API surface,
or the process rules that keep documentation current.

This page is for maintainers. If you are trying to run SpaceAGORA or find a
command, see the [User Guide](../src/index.md).

From the root of a prepared checkout, use Julia 1.12 for the maintainer check:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
```

What to read next:

- [Contributor Guide](index.md)
- [API Policy](../src/public_api_policy.md)
- [Contracts](contracts.md)

SpaceAGORA uses three documentation layers, and each one has a different owner.
Only user workflows and supported API or extension references belong in
`docs/src/`. Contributor procedures and repository contracts stay outside that
published source tree.

## 1. Generated API docs

Generated API docs are built from Julia docstrings and the explicit public API
list in:

```text
docs/public_api_symbols.jl
```

Policy:

- generated API docs are the source of truth for the stable exported package
  surface
- if root exports change, `docs/public_api_symbols.jl` must be updated in the
  same pull request
- undocumented public symbols are treated as documentation failures

## 2. Contract docs

Architecture and quality contracts remain hand-written and normative.

Sources:

- `docs/architecture/*`
- `docs/quality/*`

Policy:

- contract docs describe repository rules, ownership boundaries, and quality
  gates
- generated API docs do not replace these contract docs

## 3. Operational docs

Operational docs describe user workflows such as:

- CLI usage
- assets and licensing boundaries
- distributed/HPC launch patterns
- extension recipes and templates

Policy:

- operational docs should have either:
  - a matching test gate, or
  - a machine-readable source of truth

Examples already in the repo:

- assets are anchored to `data/assets_manifest.toml`
- public API docs are anchored to `docs/public_api_symbols.jl`
- the extension hook reference is generated from `docs/extension_hooks.jl` and
  checked by `test/gates/ci_extension_hook_reference_gate.jl`
- Wave 5 docs are anchored by `test/gates/ci_hpc_extensibility_docs_gate.jl`
- the code of the worked examples in `docs/src/user/gnc_feedback.md` and
  `docs/src/user/gnc_planning.md`, and of
  the two extension walkthroughs (`custom_effector.md`, `custom_actuator.md`) is
  the code CI runs: each block is a snippet region of a script in
  `examples/gnc_algorithm_env/` or `examples/walkthrough_env/`, compared byte for
  byte by `test/gates/ci_hpc_extensibility_docs_gate.jl`; change the page and
  the script together
- terminal commands in user docs should use cross-platform entrypoints or
  provide shell-specific variants where syntax differs

## PR review checklist

Every PR should check the docs drift items in `.github/pull_request_template.md`:

- if exports changed, update `docs/public_api_symbols.jl`
- if assets changed, update `data/assets_manifest.toml` and `docs/src/assets.md`
- if new extension points were added, update `docs/src/extensibility.md`
- if new CLI commands were added, update `docs/src/cli.md`

## Release gate

Before merge to `main`, the documentation release gate is:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
julia --project=. test/gates/ci_public_api_surface_gate.jl
julia --project=. test/gates/ci_hpc_extensibility_docs_gate.jl
```

These are treated as release blockers for the documentation surface.

## Documentation enforcement

The docs build checks the root exported public API list and fails if any symbol
listed in the generated public API page is missing documentation.

The generated undocumented-export report is written to:

```text
docs/build/undocumented_public_exports.txt
```

That report is intentionally scoped to the stable root API, not to internal
modules.

## Build the docs locally

From the repository root, use Julia 1.12 to prepare the docs environment and
build the documentation. Dependency setup may access package servers:

```text
julia --startup-file=no docs/setup.jl
julia --project=docs docs/make.jl
```

The docs environment (`docs/Project.toml`) holds Documenter and takes SpaceAGORA
from this checkout; the root project does not include Documenter. `docs/setup.jl`
generates the ignored `docs/Manifest.toml` from the root manifest on every run,
so the build uses the package versions the tests use. Run it again after the root
manifest changes.

The build checks external links, and GitHub limits anonymous requests, so a local
link check can fail with HTTP 429. Give the build a GitHub token: in bash or zsh,
`GITHUB_TOKEN=$(gh auth token) julia --project=docs docs/make.jl`; in PowerShell,
set `$env:GITHUB_TOKEN = gh auth token` first. Documenter sends the token only
with its requests to GitHub.

Open `docs/build/index.html` after a successful build. Pages under
`docs/contributing/`, `docs/architecture/` and `docs/quality/` are not rendered
into that site; review their Markdown, command context and repository links
separately when changing them.

Generated documentation stays local under the ignored `docs/build/`,
`docs/site/` and `docs/src/generated/` directories.

## Figures

Figures live in `docs/src/assets/figures/` and are committed; the documentation
build does not regenerate them. The scripts in `scripts/docs_figures/` reproduce
each one, and its README lists the page that uses each figure and what each
script needs. When a change alters what a figure shows, rerun its script and
commit the result. Diagrams are SVG, plots PNG and viewer screenshots JPEG, to
keep them small. Every image needs alt text that says what it shows.

## Regenerate common local artifacts

The repository helper recreates common ignored reports or documentation:

```text
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```

These commands execute the selected study or build; they are not required to
inspect existing simulation results. The helper's `docs` mode prepares the docs
environment and builds with it, like the direct commands above.

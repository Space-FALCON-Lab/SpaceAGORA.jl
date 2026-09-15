---
id: misc.assets_check_assets
label: check_assets
kind: function
source:
  file: src/cli/assets.jl
  symbol: check_assets
  lines:
  - 83
  - 99
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORACLI namespace supplying REPO_ROOT, the asset item constructor
    and the manifest loader.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: report
  type: AssetCheckReport
  units: n/a
  description: Repository root plus one AssetCheckItem per manifest entry, each flagged
    available/required with a licensing detail string.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
- assets
- manifest
charts:
- misc
origin: agent
---

# check_assets

## Purpose
`check_assets` answers a single operational question before a run starts: are the data roots this repository expects actually present on disk? It reads the declarative asset manifest at `data/assets_manifest.toml`, resolves every declared entry against a repository root, tests its existence, and returns an `AssetCheckReport` that the CLI prints through `render_asset_report`. Operators use it to distinguish a genuinely missing licensed dataset such as a GRAM atmosphere distribution from an optional extra that the baseline no-GRAM configuration never needs.

## Model & Assumptions
The manifest is the authority; the function invents no paths of its own. Each `[[asset]]` table declares `name`, `scope`, `relative_path`, `kind`, `required`, `licensing` and an optional `detail`. Only three kinds are modelled: `builtin` entries are unconditionally available because they ship inside the package, `directory` entries are tested with `isdir`, and `file` entries with `isfile`. Any other kind raises `ArgumentError` rather than silently reporting unavailability. A `relative_path` of `""` or `"."` resolves to the repository root itself.

## Design & Implementation
`load_asset_manifest` parses the TOML into a vector of `AssetManifestEntry` values, applying string coercion and defaults (`kind` defaults to `directory`, `required` to `false`, `licensing` to `unspecified`). `check_assets` then loops those entries, calls `_manifest_entry_path` to join `repo_root` with the split relative path, and `_manifest_entry_available` to probe the filesystem. The detail column is composed so licensing is always visible: when the entry carries no prose detail the string becomes `licensing=<value>`, otherwise the declared detail is suffixed with `[licensing=<value>]`. Both `repo_root` and `manifest_path` are keyword arguments, which lets the test suite point the check at a fixture tree.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SpaceAGORACLI namespace supplying REPO_ROOT, the asset item constructor and the manifest loader. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `report` | AssetCheckReport | n/a | — | Repository root plus one AssetCheckItem per manifest entry, each flagged available/required with a licensing detail string. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets__manifest_entry_available|_manifest_entry_available]] · `callees` → `callers` · feedback · `src/cli/assets.jl:79-79`
- [[cli.assets_assetcheckitem|AssetCheckItem]] · `callees` → `callers` · call · `src/cli/assets.jl:20-20`
- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · feedback · `src/cli/assets.jl:136-136`
- [[cli.run_cli|run_cli]] · `callees` → `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
- [[spaceagora.spaceagora_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:547-547`

**Downstream**

- `callees` → [[cli.assets__asset_item|_asset_item]] · `callers` · feedback · `src/cli/assets.jl:88-88`
- `callees` → [[cli.assets__manifest_entry_available|_manifest_entry_available]] · `callers` · call · `src/cli/assets.jl:92-92`
- `callees` → [[cli.assets__manifest_entry_path|_manifest_entry_path]] · `callers` · call · `src/cli/assets.jl:91-91`
- `callees` → [[cli.assets_assetcheckreport|AssetCheckReport]] · `callers` · call · `src/cli/assets.jl:98-98`
- `callees` → [[cli.assets_load_asset_manifest|load_asset_manifest]] · `callers` · call · `src/cli/assets.jl:84-84`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/cli/assets.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
Availability means existence, not usability: the check never opens a file, verifies a checksum, inspects permissions, or confirms that a directory contains the expected contents. An empty directory of the right name reports as available. The report is a snapshot with no timestamp, so a stale printout can outlive the filesystem state it described. A malformed or absent manifest raises during parsing instead of producing a partial report.

## Provenance
Mapped from `src/cli/assets.jl:83-99`, with supporting definitions at lines 8-76.

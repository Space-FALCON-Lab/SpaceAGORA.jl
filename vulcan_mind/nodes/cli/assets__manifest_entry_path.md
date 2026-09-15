---
id: cli.assets__manifest_entry_path
label: _manifest_entry_path
kind: function
source:
  file: src/cli/assets.jl
  symbol: _manifest_entry_path
  lines:
  - 65
  - 65
inputs:
- id: repo_root
  type: String
  units: n/a
  required: true
  description: Positional argument `repo_root`.
- id: entry
  type: AssetManifestEntry
  units: n/a
  required: true
  description: Positional argument `entry`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Any
  units: n/a
  description: Return value of `_manifest_entry_path`. Returns `joinpath(repo_root,
    splitpath(entry.relative_path)...)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# _manifest_entry_path

## Purpose
Resolves the absolute filesystem location of a manifest entry by combining the repository root with the entry's declared `relative_path`. It is the single place where the special-case meaning of an empty or `.` path (the repo root itself) is implemented.

## Design & Implementation
`@inline` function taking `repo_root::String` and `entry::AssetManifestEntry`. If `entry.relative_path` is `""` or `"."` it returns `repo_root` unchanged. Otherwise it returns `joinpath(repo_root, splitpath(entry.relative_path)...)`, splitting the manifest path into components first so that a manifest written with forward slashes is reassembled with the host separator. It is called by `_manifest_entry_available` before probing the filesystem and by `check_assets` to populate `AssetCheckItem.path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `repo_root` | String | n/a | yes | Positional argument `repo_root`. |
| in | `entry` | AssetManifestEntry | n/a | yes | Positional argument `entry`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_manifest_entry_path`. Returns `joinpath(repo_root, splitpath(entry.relative_path)...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets__manifest_entry_available|_manifest_entry_available]] · `callees` → `callers` · call · `src/cli/assets.jl:71-71`
- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:91-91`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No normalisation such as `normpath` or `abspath` is applied, so a `relative_path` containing `..` or an absolute path is joined verbatim and can point outside the repository. Trailing separators or mixed separators survive `splitpath` in host-specific ways. The function does not check that `repo_root` exists.

## Provenance
Mapped from `src/cli/assets.jl` line 65.

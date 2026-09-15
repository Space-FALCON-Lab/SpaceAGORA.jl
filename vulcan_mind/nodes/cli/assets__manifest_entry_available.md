---
id: cli.assets__manifest_entry_available
label: _manifest_entry_available
kind: function
source:
  file: src/cli/assets.jl
  symbol: _manifest_entry_available
  lines:
  - 70
  - 70
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
  description: Return value of `_manifest_entry_available`.
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

# _manifest_entry_available

## Purpose
Determines whether a manifest entry's asset is present on disk, choosing the probe by the entry's declared `kind`. It supplies the `available` flag for every `AssetCheckItem` produced by `check_assets`.

## Design & Implementation
`@inline` function of `repo_root::String` and `entry::AssetManifestEntry`. It first computes `path = _manifest_entry_path(repo_root, entry)`, then dispatches on the string `entry.kind`: `"builtin"` returns `true` unconditionally (the asset ships with the package), `"directory"` returns `isdir(path)`, and `"file"` returns `isfile(path)`. Any other kind falls through to `throw(ArgumentError("Unsupported asset manifest kind '<kind>' for <name>."))`, naming both the offending kind and the entry.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `repo_root` | String | n/a | yes | Positional argument `repo_root`. |
| in | `entry` | AssetManifestEntry | n/a | yes | Positional argument `entry`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_manifest_entry_available`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:92-92`

**Downstream**

- `callees` → [[cli.assets__manifest_entry_path|_manifest_entry_path]] · `callers` · call · `src/cli/assets.jl:71-71`
- `callees` → [[misc.assets_check_assets|check_assets]] · `callers` · feedback · `src/cli/assets.jl:79-79`
- `callees` → [[spaceagora.spaceagora_check_assets|check_assets]] · `callers` · call · `src/cli/assets.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
The path is resolved even for `"builtin"` entries, which is harmless but wasted work. A directory that exists but is empty or lacks its expected contents counts as available; no manifest field allows specifying a sentinel file inside the directory. Permission-denied conditions make `isdir`/`isfile` return `false`, so an unreadable asset is reported as missing rather than inaccessible. Because the kind check is a chain of string comparisons, a capitalised kind such as `"Directory"` throws.

## Provenance
Mapped from `src/cli/assets.jl` line 70.

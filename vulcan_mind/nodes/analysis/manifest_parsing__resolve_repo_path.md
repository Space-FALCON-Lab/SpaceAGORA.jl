---
id: analysis.manifest_parsing__resolve_repo_path
label: _resolve_repo_path
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _resolve_repo_path
  lines:
  - 426
  - 426
inputs:
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
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
  type: String
  units: n/a
  description: Return value of `_resolve_repo_path`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _resolve_repo_path

## Purpose
Anchors a relative manifest path at the repository root so scenarios are runnable from any working directory.

## Design & Implementation
Returns the path unchanged if `isabspath`, otherwise `normpath(joinpath(REPO_ROOT, path))`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_resolve_repo_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:548-548`
- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:378-378`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Paths are relative to the repository root, not to the manifest file, so a manifest stored outside the repo cannot use relative references to its own neighbours.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 426.

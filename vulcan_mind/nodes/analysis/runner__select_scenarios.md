---
id: analysis.runner__select_scenarios
label: _select_scenarios
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _select_scenarios
  lines:
  - 341
  - 341
inputs:
- id: scenarios
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `scenarios`.
- id: requested
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `requested`.
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
  type: AbstractArray
  units: n/a
  description: Return value of `_select_scenarios`. Returns `[sc for sc in scenarios
    if lowercase(String(sc.name)) in wanted]`.
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

# _select_scenarios

## Purpose
`_select_scenarios` filters the manifest's scenario list down to the names requested on the command line or in a `VerificationRequest`, failing loudly on unknown names so that a typo cannot silently skip a CI scenario. An empty request keeps every scenario.

## Design & Implementation
Signature `(scenarios::AbstractVector, requested::Vector{String})`. It normalises `requested` by `strip` and `lowercase`, dropping empty strings and deduplicating with `unique`, and returns `scenarios` unchanged if nothing remains. Manifest names are lowercased with `lowercase(String(sc.name))`; `setdiff(wanted, names)` yields unknown requests, and a non-empty set throws `ArgumentError` listing both the unknown names and all manifest names. The result is a comprehension preserving manifest order for scenarios whose lowercased name is in `wanted`. Nothing is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scenarios` | AbstractVector | n/a | yes | Positional argument `scenarios`. |
| in | `requested` | Vector{String} | n/a | yes | Positional argument `requested`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `_select_scenarios`. Returns `[sc for sc in scenarios if lowercase(String(sc.name)) in wanted]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__final_run_or_reused_eval|_final_run_or_reused_eval]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:335-335`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:355-355`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Matching is exact after case-folding and trimming; there is no glob or prefix matching. The returned order is manifest order, not request order, so callers cannot control execution sequence by argument order. If the manifest itself contains two scenarios differing only by case both are selected. `sc.name` must be convertible with `String(...)`, otherwise a `MethodError` surfaces here rather than in manifest loading.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 341.

---
id: gnc.target_energy_bracketing__edg_validate_symbol_set
label: _edg_validate_symbol_set
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_validate_symbol_set
  lines:
  - 10
  - 10
inputs:
- id: values
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `values`.
- id: allowed
  type: Set{Symbol}
  units: n/a
  required: true
  description: Positional argument `allowed`.
- id: label
  type: String
  units: n/a
  required: true
  description: Positional argument `label`.
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
  description: Return value of `_edg_validate_symbol_set`. Returns `values`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _edg_validate_symbol_set

## Purpose
`_edg_validate_symbol_set` checks that a tuple of mode symbols is non-empty and that every member belongs to an allowed set, raising a descriptive `ArgumentError` otherwise. The `AerobrakingEnergyDepletionConfig` constructor uses it for `guidance_modes` (against `_EDG_GUIDANCE_MODES`) and `max_energy_submodes` (against `_EDG_MAX_ENERGY_SUBMODES`).

## Design & Implementation
Declared `@inline` with signature `(values::Tuple, allowed::Set{Symbol}, label::String)`. It first throws `ArgumentError("$(label) must not be empty.")` when `values` is empty; then for each `value` it tests `value in allowed` and throws `ArgumentError` naming the label, the offending value, and `sort!(collect(allowed))` so the message lists every permitted symbol in order. On success it returns `values` unchanged, enabling the constructor to bind the validated tuple directly. Nothing is mutated except the temporary sorted copy built for the error message.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Tuple | n/a | yes | Positional argument `values`. |
| in | `allowed` | Set{Symbol} | n/a | yes | Positional argument `allowed`. |
| in | `label` | String | n/a | yes | Positional argument `label`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_validate_symbol_set`. Returns `values`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.target_energy_bracketing_aerobrakingenergydepletionconfig|AerobrakingEnergyDepletionConfig]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:63-63`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The allowed sets are module constants (`:max_energy_depletion`, `:targeting`; `:heat_rate`, `:structural_load`, `:heat_load`), so extending the mode vocabulary requires editing this file. Validation is per-element only: combinations that are semantically inconsistent (for example `:targeting` without any submodes that can bracket energy) are not detected. Duplicates within `values` are accepted. The function does not verify element types beyond membership, relying on `_edg_symbol_tuple` having produced symbols.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 10.

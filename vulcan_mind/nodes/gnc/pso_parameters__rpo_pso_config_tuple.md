---
id: gnc.pso_parameters__rpo_pso_config_tuple
label: _rpo_pso_config_tuple
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: _rpo_pso_config_tuple
  lines:
  - 310
  - 310
inputs:
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: NamedTuple
  units: n/a
  description: Return value of `_rpo_pso_config_tuple`. Returns `NamedTuple{names}(map(name
    -> getfield(cfg, name), names))`.
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

# _rpo_pso_config_tuple

## Purpose
Converts an `RPOPSOConfig` into a `NamedTuple` with identical field names and values so it can be merged with keyword overrides and splatted back into the keyword constructor.

## Design & Implementation
Uses `fieldnames(RPOPSOConfig)` to obtain the field-name tuple, maps `getfield(cfg, name)` over it, and wraps the result as `NamedTuple{names}(...)`. Because `fieldnames` is evaluated at runtime on the type, the function automatically tracks any change to the struct's field list. It is called by `rpo_pso_config` as the base of the `merge` chain.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NamedTuple | n/a | — | Return value of `_rpo_pso_config_tuple`. Returns `NamedTuple{names}(map(name -> getfield(cfg, name), names))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:584-584`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Allocates a tuple of roughly 115 elements per call, and the `map` over a runtime `names` tuple may not be fully type-inferred, so this is not intended for hot loops. Field order is preserved but nothing prevents a caller from constructing a NamedTuple with extra keys that `RPOPSOConfig(; values...)` will then reject with a `MethodError`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 310.

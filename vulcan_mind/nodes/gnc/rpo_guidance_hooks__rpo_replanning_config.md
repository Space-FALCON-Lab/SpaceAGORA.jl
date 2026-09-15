---
id: gnc.rpo_guidance_hooks__rpo_replanning_config
label: _rpo_replanning_config
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: _rpo_replanning_config
  lines:
  - 60
  - 60
inputs:
- id: model
  type: RPOGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: RPOReplanningConfig
  units: n/a
  description: Return value of `_rpo_replanning_config`. Returns `RPOReplanningConfig(;
    model.replanning_config...)`.
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

# _rpo_replanning_config

## Purpose
Normalises the replanning configuration attached to an RPO guidance model into a concrete `RPOReplanningConfig`, so callers can accept either a fully built configuration object or a loose keyword collection from a scenario file.

## Design & Implementation
Returns `nothing` when `model.replanning_config === nothing`, signalling that replanning is entirely absent. If the stored value already `isa RPOReplanningConfig` it is returned unchanged. Otherwise the value is splatted as keyword arguments into the constructor, `RPOReplanningConfig(; model.replanning_config...)`, which accepts a named tuple or a dictionary of symbol keys and fills the remaining fields from their defaults.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOReplanningConfig | n/a | — | Return value of `_rpo_replanning_config`. Returns `RPOReplanningConfig(; model.replanning_config...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:84-84`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`

**Downstream**

- `callees` → [[gnc.replanning_rporeplanningconfig|RPOReplanningConfig]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations
The splat path reconstructs the configuration object on every call rather than caching it, and `maybe_update_rpo_replanning!` calls this once per guidance step, so a scenario that stores a named tuple pays a constructor call every step. An unrecognised key in the stored collection raises a `MethodError` from the keyword constructor rather than a message naming the offending field, and a non-splattable stored value fails the same way. The distinction between a missing configuration and a present-but-disabled one is left to the caller, which checks `config.enabled` separately.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 60.

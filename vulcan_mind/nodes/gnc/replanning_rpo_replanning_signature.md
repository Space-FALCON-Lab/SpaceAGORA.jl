---
id: gnc.replanning_rpo_replanning_signature
label: rpo_replanning_signature
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_replanning_signature
  lines:
  - 127
  - 127
inputs:
- id: spheres
  type: AbstractVector{RPOReplanningSphere}
  units: n/a
  required: true
  description: Positional argument `spheres`.
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
  description: Return value of `rpo_replanning_signature`. Returns `acc`.
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

# rpo_replanning_signature

## Purpose
Reduces the active obstacle set to a single `UInt` fingerprint so the guidance loop can detect whether the obstacle configuration has materially changed since the last plan without comparing sphere lists element by element.

## Design & Implementation
Signature `rpo_replanning_signature(spheres::AbstractVector{RPOReplanningSphere})`. For an empty vector it returns `UInt(0)`. Otherwise it starts from the constant `0x9e3779b97f4a7c15` (the 64-bit golden-ratio constant) and XORs in `hash((label, round.(Tuple(center_rtn); digits=4), round(radius_m; digits=4)))` for each sphere. Rounding centre and radius to four decimals (0.1 mm) suppresses spurious changes from floating-point noise in drifted centres. The XOR fold makes the result independent of sphere ordering.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spheres` | AbstractVector{RPOReplanningSphere} | n/a | yes | Positional argument `spheres`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_replanning_signature`. Returns `acc`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:94-94`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
XOR folding means two identical spheres cancel to the seed value, so a set containing a duplicated obstacle hashes the same as the set without either copy. Lifetime and velocity fields are excluded, so two spheres with identical current centre but different drift produce the same signature. `hash` of tuples is not stable across Julia versions, so signatures must not be persisted. Because drifting spheres change centre every step, a moving obstacle changes the signature continuously and defeats the intent of a stable key.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 127.

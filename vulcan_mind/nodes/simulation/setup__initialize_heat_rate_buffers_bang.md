---
id: simulation.setup__initialize_heat_rate_buffers_bang
label: _initialize_heat_rate_buffers!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_heat_rate_buffers!
  lines:
  - 127
  - 127
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Nothing
  units: n/a
  description: Return value of `_initialize_heat_rate_buffers!`; mutates `p` in place.
    Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _initialize_heat_rate_buffers!

## Purpose
Sizes and zeros the per-spacecraft, per-link heat-rate scratch vectors in `p.shared_buffers.heat_rates` so thermal callbacks can write into preallocated storage during integration.

## Design & Implementation
Computes `n_sats = length(p.args.dynamics_model.spacecraft)` and `resize!`s the outer vector if its length differs. For each satellite `i` it takes `n_links = length(spacecraft[i].links)`; if the existing slot is not a `Vector{Float64}` (for example `#undef` after a resize) it allocates `Float64[]` and stores it, then `resize!`s that inner vector to `n_links` and `fill!`s it with `0.0`. The loop is `@inbounds`. Returns `nothing`. Mutates `p.shared_buffers.heat_rates` and its elements in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_heat_rate_buffers!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:190-190`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`rates isa Vector{Float64}` on an `#undef` slot throws `UndefRefError` before the check can replace it; the code only works because `resize!` growth leaves slots undefined only if the eltype is not isbits, so the exact behaviour depends on how `heat_rates` was declared. Shrinking `n_sats` leaves stale inner vectors unreachable but not freed. Not thread-safe against concurrent readers of the buffers.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 127.

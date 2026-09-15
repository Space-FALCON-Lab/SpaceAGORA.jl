---
id: simulation.effector_sampling__buffered_atmosphere_valid
label: _buffered_atmosphere_valid
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _buffered_atmosphere_valid
  lines:
  - 127
  - 127
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Bool
  units: n/a
  description: Return value of `_buffered_atmosphere_valid`.
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

# _buffered_atmosphere_valid

## Purpose
Decides whether the shared atmosphere buffer for a satellite may be reused at time `t`, under either exact-time or freeze-per-step semantics.

## Design & Implementation
Returns false if the satellite index exceeds the buffer length. In freeze-per-step mode it requires only that the buffered sample time be finite, trusting the density callback's once-per-accepted-step resample for every RHS stage. Otherwise it requires the buffered time to equal `t` exactly. `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_buffered_atmosphere_valid`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling_sample_buffered_atmosphere|sample_buffered_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:142-142`

**Downstream**

- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:135-135`
<!-- vulcan:connections:end -->

## Limitations
Exact floating-point equality on `t` means a stage evaluated at a time computed by a slightly different arithmetic path misses the buffer and pays a fresh sample; this is deliberate conservatism but can silently double density cost.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 127.

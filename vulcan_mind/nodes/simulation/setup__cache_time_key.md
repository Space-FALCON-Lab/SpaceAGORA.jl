---
id: simulation.setup__cache_time_key
label: _cache_time_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _cache_time_key
  lines:
  - 240
  - 240
inputs:
- id: x
  type: Float64
  units: n/a
  required: true
  description: Positional argument `x`.
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
  type: Int64
  units: n/a
  description: Return value of `_cache_time_key`.
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

# _cache_time_key

## Purpose
Quantises a floating-point time in seconds to an `Int64` microsecond count so that ephemeris reuse keys compare exactly instead of suffering from floating-point representation noise.

## Theory & Math
The key is $k = \operatorname{round}(10^{6} x)$ where $x$ is a time in seconds and $k$ is the resulting integer microsecond count.

## Design & Implementation
`_cache_time_key(x::Float64)::Int64` returns `round(Int64, x * 1e6)`. Applied to `et_start`, `mission_end_s`, and `dt_s` in every `_*_ephemeris_reuse_key` builder. Rounding (not truncation) means values within half a microsecond of each other collide, which is the intended equivalence.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Float64 | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int64 | n/a | — | Return value of `_cache_time_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:283-283`
- [[simulation.setup__planet_frame_ephemeris_reuse_key|_planet_frame_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:294-294`
- [[simulation.setup__srp_ephemeris_reuse_key|_srp_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:273-273`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Ephemeris times around J2000 in seconds (order 1e9) times 1e6 stay well inside `Int64`, but inputs beyond about 9.2e12 s overflow and `round` throws `InexactError`. `NaN` or `Inf` also throw. Two runs differing by less than 1 μs in `dt_s` share caches even if intended to differ.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 240.

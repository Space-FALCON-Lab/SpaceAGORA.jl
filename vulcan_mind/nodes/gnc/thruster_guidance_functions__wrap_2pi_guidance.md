---
id: gnc.thruster_guidance_functions__wrap_2pi_guidance
label: _wrap_2pi_guidance
kind: function
source:
  file: src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl
  symbol: _wrap_2pi_guidance
  lines:
  - 14
  - 14
inputs:
- id: theta
  type: Float64
  units: n/a
  required: true
  description: Positional argument `θ`.
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
  type: Float64
  units: n/a
  description: Return value of `_wrap_2pi_guidance`.
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

# _wrap_2pi_guidance

## Purpose
Normalises an angle in radians into the half-open interval [0, 2π) so true-anomaly comparisons in the apoapsis-targeting guidance can test pre/post-apoapsis with a simple `ν <= π` check.

## Design & Implementation
An `@inline` function `_wrap_2pi_guidance(θ::Float64)::Float64`. It computes `θw = mod(θ, 2pi)`, which in Julia already returns a value with the sign of the divisor (non-negative), and then defensively adds `2pi` if `θw < 0.0`. The result is used in `calcGuidanceEffect!` for `ApoapsisTargetPeriapsisRaiseGuidanceModel` to compute `distance_to_apoapsis = π - ν`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `theta` | Float64 | n/a | yes | Positional argument `θ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_wrap_2pi_guidance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:173-173`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The negative-branch guard is unreachable for finite inputs because `mod` with a positive divisor is already non-negative, so it is dead code. `mod(θ, 2pi)` on a NaN or Inf input returns NaN, which then fails every downstream comparison and silently suppresses the burn command. Floating-point `mod` near multiples of 2π can return a value equal to `2pi` minus an ulp rather than 0.

## Provenance
Mapped from `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` line 14.

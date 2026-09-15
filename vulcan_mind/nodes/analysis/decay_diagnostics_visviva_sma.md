---
id: analysis.decay_diagnostics_visviva_sma
label: visviva_sma
kind: function
source:
  file: src/analysis/verification/telemetry_verification/decay_diagnostics.jl
  symbol: visviva_sma
  lines:
  - 22
  - 22
inputs:
- id: r_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `r_m`.
- id: v_mps
  type: Real
  units: n/a
  required: true
  description: Positional argument `v_mps`.
- id: mu
  type: Real
  units: n/a
  required: true
  description: Positional argument `mu`.
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
  description: Return value of `visviva_sma`. Returns `1.0 / (2.0 / r_m - v_mps^2
    / mu)`.
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

# visviva_sma

## Purpose

`visviva_sma(r_m, v_mps, mu)` returns the osculating semi-major axis in metres from the instantaneous position magnitude `r_m` in metres, the speed `v_mps` in metres per second and the gravitational parameter `mu`. It converts a telemetry or simulation state series into the SMA series that the decay estimators regress.

## Design & Implementation

There are two methods. The scalar method is marked `@inline` and evaluates `1.0 / (2.0 / r_m - v_mps^2 / mu)` directly, so it allocates nothing and inlines into hot loops. The vector method takes `AbstractVector` arguments for `r_m` and `v_mps` with a scalar `mu` and broadcasts the scalar method elementwise, returning a new vector of the same length. Both accept any `Real` inputs and promote through the literal `Float64` constants in the expression.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_m` | Real | n/a | yes | Positional argument `r_m`. |
| in | `v_mps` | Real | n/a | yes | Positional argument `v_mps`. |
| in | `mu` | Real | n/a | yes | Positional argument `mu`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `visviva_sma`. Returns `1.0 / (2.0 / r_m - v_mps^2 / mu)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`

**Downstream**

- `callees` → [[envana.ana_decay_diagnostics_secular_sma_slope|secular_sma_slope]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:27-27`
<!-- vulcan:connections:end -->

## Theory & Math
This is the vis-viva relation solved for the semi-major axis:

$$v^2 = \mu\left(\frac{2}{r} - \frac{1}{a}\right) \quad\Longrightarrow\quad a = \left(\frac{2}{r} - \frac{v^2}{\mu}\right)^{-1}.$$

## Limitations

For a parabolic state the bracket goes to zero and the result diverges; for hyperbolic states it returns a negative semi-major axis rather than raising. `r_m` must be non-zero, and the two vectors must have matching lengths, which the broadcast enforces only through its own shape error. Because the value is osculating, it carries the full short-period J2 signature and must not be regressed for a secular slope without the zero-reference correction described in the module header.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/decay_diagnostics.jl` line 22.

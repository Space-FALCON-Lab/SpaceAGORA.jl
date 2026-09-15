---
id: analysis.comparison_metrics__rates
label: _rates
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _rates
  lines:
  - 55
  - 55
inputs:
- id: axis
  type: Any
  units: n/a
  required: true
  description: Positional argument `axis`.
- id: val
  type: Any
  units: n/a
  required: true
  description: Positional argument `val`.
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
  description: Return value of `_rates`. Returns `(         (val[2:end] .- val[1:(end
    - 1)]) ./ (axis[2:end] .- axis[1:(end - 1)])`.
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

# _rates

## Purpose

A closure defined inside `_apo_decay_diagnostic` that converts a sampled altitude series into finite-difference decay rates together with the interval midpoints at which those rates apply. It is the shared step that puts telemetry apoapsis and simulated apoapsis onto comparable per-orbit rate curves.

## Design & Implementation

`_rates(axis, val)` returns a two-element tuple. The first element is the backward difference quotient `(val[2:end] .- val[1:end-1]) ./ (axis[2:end] .- axis[1:end-1])`, giving altitude change per unit of axis (kilometres per orbit for the apoapsis rows). The second is the midpoint axis `(axis[2:end] .+ axis[1:end-1]) ./ 2.0`, one value shorter than the input, which the caller uses as the common abscissa when interpolating the simulated rate onto the telemetry intervals.

## Theory & Math

For samples $(a_k, v_k)$ the closure computes, for each interval $k$,

$$\dot{v}_k = \frac{v_{k+1} - v_k}{a_{k+1} - a_k}, \qquad \bar{a}_k = \frac{a_{k+1} + a_k}{2}$$

where $a_k$ is the orbit index and $v_k$ the apoapsis altitude in km, so $\dot{v}_k$ carries units of km per orbit and $\bar{a}_k$ is the interval-centred abscissa. Centring at the midpoint makes the difference quotient second-order accurate for a smooth decay profile.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `axis` | Any | n/a | yes | Positional argument `axis`. |
| in | `val` | Any | n/a | yes | Positional argument `val`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rates`. Returns `(         (val[2:end] .- val[1:(end - 1)]) ./ (axis[2:end] .- axis[1:(end - 1)])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl`

**Downstream**

- `callees` → [[analysis.comparison_metrics__interp_linear|_interp_linear]] · `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:77-77`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations

No guard exists against a zero axis increment, so repeated orbit indices produce infinite or `NaN` rates that then flow into the ratio statistics. Each call allocates four temporary arrays through the slicing and broadcast, which is acceptable at the small orbit counts involved but not for dense series. Being a local closure, it cannot be reused or tested outside `_apo_decay_diagnostic`, and it inherits that function's assumption that the axis is monotonically increasing.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/comparison_metrics.jl` line 55.

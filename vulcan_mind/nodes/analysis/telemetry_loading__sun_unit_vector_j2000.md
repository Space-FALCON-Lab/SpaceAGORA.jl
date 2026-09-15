---
id: analysis.telemetry_loading__sun_unit_vector_j2000
label: _sun_unit_vector_j2000
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _sun_unit_vector_j2000
  lines:
  - 398
  - 398
inputs:
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: el_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_s`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_sun_unit_vector_j2000`.
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

# _sun_unit_vector_j2000

## Purpose
Cheap analytic Sun direction in J2000 coordinates at a given elapsed time past the scenario epoch, used only to screen telemetry samples by dayside or nightside illumination where roughly 0.01 degree accuracy is sufficient and a SPICE ephemeris call per sample would be unnecessary.

## Theory & Math
With $n$ = days past J2000.0:
$$L = 280.460^\circ + 0.9856474^\circ n,\quad g = 357.528^\circ + 0.9856003^\circ n$$
$$\lambda = L + 1.915^\circ \sin g + 0.020^\circ \sin 2g,\quad \varepsilon = 23.439^\circ - 4\times10^{-7\,\circ} n$$
$$\hat s = (\cos\lambda,\ \cos\varepsilon \sin\lambda,\ \sin\varepsilon \sin\lambda)$$

## Design & Implementation
Builds a `DateTime` from `initial_time` fields, adding the fractional seconds as `Millisecond(round(Int, 1000*second))`, then computes days since J2000 as `n_days = (epoch - 2000-01-01T12:00) / 86400e3 + el_s / 86400`. It evaluates the Meeus low-precision solar model: mean longitude `L = 280.460 + 0.9856474 n` deg, mean anomaly `g = 357.528 + 0.9856003 n` deg, ecliptic longitude `λ = L + 1.915 sin g + 0.020 sin 2g` deg, obliquity `ε = 23.439 - 4e-7 n` deg, and returns the unit vector `(cos λ, cos ε sin λ, sin ε sin λ)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `el_s` | Float64 | n/a | yes | Positional argument `el_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_sun_unit_vector_j2000`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:338-338`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:402-402`
<!-- vulcan:connections:end -->

## Limitations
The model is geocentric: it gives the Earth-to-Sun direction and is only approximately correct as a Sun direction for spacecraft at other planets, so the mask is questionable for non-Earth scenarios. `mod(..., 360)` is applied to `L` and `g` but not to `λ`, which is harmless for trig but means the intermediate is unbounded. Fractional seconds are rounded to milliseconds. Accuracy degrades outside roughly 1950 to 2050 as the linear terms drift.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 398.

---
id: envana.ana_telemetry_loading_transform_state
label: _transform_state
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _transform_state
  lines:
  - 35
  - 39
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace providing the SPICE sxform binding
    and StaticArrays types.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: state
  type: Tuple{SVector{3,Float64},SVector{3,Float64}}
  units: m,m/s
  description: Position and velocity vectors rotated into the destination reference
    frame.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# _transform_state

## Purpose
`_transform_state` converts a full position and velocity pair between two SPICE-named reference frames at a given epoch, which is the operation telemetry loading needs whenever recorded data arrives in a planet-fixed frame but the comparison is performed in an inertial one.

## Theory & Math
A six-by-six state transformation matrix is required rather than a three-by-three rotation because velocity in a rotating frame picks up the transport term. Writing the block form as `[[R, 0], [dR/dt, R]]`, the transformed state is `r' = R r` and `v' = (dR/dt) r + R v`, which is exactly the Coriolis relation `v' = R (v + omega x r)` when `R` is a rotation with angular velocity `omega` in rad/s. Position enters in metres and velocity in metres per second, and both units are preserved because the transformation is a pure rotation with a derivative block of units per second.

## Model & Assumptions
Frame names are passed as strings and resolved by SPICE, so both must be known to the loaded kernel pool; `et` is ephemeris time in TDB seconds past J2000. The function assumes the caller has already furnished the necessary frame kernels and, because CSPICE is not reentrant, that the surrounding code holds the shared SPICE lock. No unit conversion is applied, so the caller must supply metres and metres per second rather than the kilometre-based units SPICE itself commonly uses.

## Design & Implementation
The body is four lines. `sxform(from_frame, to_frame, et)` returns the transformation, which is immediately wrapped as an `SMatrix{6,6,Float64}` so the subsequent multiply is a statically sized operation with no heap allocation. The six state components are packed into an `SVector{6,Float64}`, multiplied, and then unpacked into two `SVector{3,Float64}` results, keeping the whole call stack-resident. The sibling `_planet_fixed_to_j2000_state` immediately below wraps this with the planet-fixed frame name lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace providing the SPICE sxform binding and StaticArrays types. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `state` | Tuple{SVector{3,Float64},SVector{3,Float64}} | m,m/s | — | Position and velocity vectors rotated into the destination reference frame. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.telemetry_loading__j2000_to_planet_fixed_state|_j2000_to_planet_fixed_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:56-56`
- [[analysis.telemetry_loading__planet_fixed_to_j2000_state|_planet_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:45-45`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Every call crosses into native SPICE, so per-sample use over a long telemetry file is measurably slower than caching the matrix for repeated epochs. No validation is performed on the frame names, so a typo surfaces as a SPICE error rather than a Julia-level message. The function is not thread-safe on its own; the lock must be held by the caller.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/telemetry_loading.jl:35-39`, together with its immediate caller `_planet_fixed_to_j2000_state` at line 41 of the same file.

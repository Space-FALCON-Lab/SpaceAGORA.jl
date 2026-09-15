---
id: analysis.example_support_make_three_body_spacecraft
label: make_three_body_spacecraft
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: make_three_body_spacecraft
  lines:
  - 61
  - 61
inputs:
- id: bus_dims
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Keyword argument `bus_dims`.
- id: panel_dims
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Keyword argument `panel_dims`.
- id: bus_mass
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `bus_mass`.
- id: panel_mass_each
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `panel_mass_each`.
- id: panel_offset_y
  type: Float64
  units: n/a
  required: true
  description: Keyword argument `panel_offset_y`.
- id: ic
  type: SM.AbstractInitialCondition
  units: n/a
  required: true
  description: Keyword argument `ic`.
- id: reflection_coefficient
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `reflection_coefficient` (default `1.0`).
- id: prop_mass
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `prop_mass` (default `0.0`).
- id: id
  type: Int64
  units: n/a
  required: false
  description: Keyword argument `id` (default `1`).
- id: bus_ram_face
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `bus_ram_face` (default `:legacy`).
- id: bus_attitude_q
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Keyword argument `bus_attitude_q` (default `nothing`).
- id: panel_attitude_q_left
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Keyword argument `panel_attitude_q_left` (default `nothing`).
- id: panel_attitude_q_right
  type: Union{Nothing, NTuple{4, Float64}}
  units: n/a
  required: false
  description: Keyword argument `panel_attitude_q_right` (default `nothing`).
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
  type: SM.SpacecraftModel
  units: n/a
  description: Return value of `make_three_body_spacecraft`. Returns `SM.SpacecraftModel(`.
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

# make_three_body_spacecraft

## Purpose
Builds the standard three-link example vehicle — a bus with two solar panel wings — with the reference-area conventions the free-molecular aerodynamics expects.

## Design & Implementation
Validates `bus_ram_face` as `:legacy` or `:frontal` and computes the bus reference area as `dims[1]*dims[3]` or `dims[2]*dims[3]` respectively; `:legacy` preserves the historical value so every previously calibrated scenario is bit-for-bit unchanged, while `:frontal` matches the face normal to the flow that the Hart coefficients are normalised by. Two panel links are created, offset by `±panel_offset_y` along y, each carrying `panel_dims[2]*panel_dims[3]` as its area — so `panel_dims[2]` must be the per-wing half-span. Attitudes come through `_link_q`. The `SpacecraftModel` is assembled with no joints, the bus as root, total mass as the sum of the three links, and the caller's `prop_mass`, initial condition and `id`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bus_dims` | NTuple{3, Float64} | n/a | yes | Keyword argument `bus_dims`. |
| in | `panel_dims` | NTuple{3, Float64} | n/a | yes | Keyword argument `panel_dims`. |
| in | `bus_mass` | Float64 | n/a | yes | Keyword argument `bus_mass`. |
| in | `panel_mass_each` | Float64 | n/a | yes | Keyword argument `panel_mass_each`. |
| in | `panel_offset_y` | Float64 | n/a | yes | Keyword argument `panel_offset_y`. |
| in | `ic` | SM.AbstractInitialCondition | n/a | yes | Keyword argument `ic`. |
| in | `reflection_coefficient` | Float64 | n/a | no | Keyword argument `reflection_coefficient` (default `1.0`). |
| in | `prop_mass` | Float64 | n/a | no | Keyword argument `prop_mass` (default `0.0`). |
| in | `id` | Int64 | n/a | no | Keyword argument `id` (default `1`). |
| in | `bus_ram_face` | Symbol | n/a | no | Keyword argument `bus_ram_face` (default `:legacy`). |
| in | `bus_attitude_q` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Keyword argument `bus_attitude_q` (default `nothing`). |
| in | `panel_attitude_q_left` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Keyword argument `panel_attitude_q_left` (default `nothing`). |
| in | `panel_attitude_q_right` | Union{Nothing, NTuple{4, Float64}} | n/a | no | Keyword argument `panel_attitude_q_right` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SM.SpacecraftModel | n/a | — | Return value of `make_three_body_spacecraft`. Returns `SM.SpacecraftModel(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_spacecraft|_make_spacecraft]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:351-351`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:4-4`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Passing the full array span as `panel_dims[2]` silently doubles the drag and SRP area — the exact defect the April 2026 examples carried until August 2026 — and nothing here can detect it. The inertia passed to the model is the bus's alone, so the panels contribute mass but no rotational inertia.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 61.

---
id: gnc.planner_comparison_rpo_lqmpc_tracking_fuel_used_pct
label: rpo_lqmpc_tracking_fuel_used_pct
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_lqmpc_tracking_fuel_used_pct
  lines:
  - 456
  - 456
inputs:
- id: fuel_used_kg
  type: Real
  units: n/a
  required: true
  description: Positional argument `fuel_used_kg`.
- id: tracking
  type: RPOLQMPCTrackingSettings
  units: n/a
  required: true
  description: Positional argument `tracking`.
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
  description: Return value of `rpo_lqmpc_tracking_fuel_used_pct`.
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

# rpo_lqmpc_tracking_fuel_used_pct

## Purpose
Converts consumed propellant mass in kilograms to a percentage of the configured `propellant_mass_kg`, the fuel metric reported in comparison CSVs and plots.

## Design & Implementation
`rpo_lqmpc_tracking_fuel_used_pct(fuel_used_kg::Real, tracking::RPOLQMPCTrackingSettings)::Float64` reads `propellant_mass = Float64(tracking.propellant_mass_kg)` and returns `100 * fuel_used_kg / propellant_mass` when that mass is finite and positive, otherwise `NaN`. Pure and allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `fuel_used_kg` | Real | n/a | yes | Positional argument `fuel_used_kg`. |
| in | `tracking` | RPOLQMPCTrackingSettings | n/a | yes | Positional argument `tracking`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `rpo_lqmpc_tracking_fuel_used_pct`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_track_retimed_path_lqmpc|rpo_track_retimed_path_lqmpc]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:526-526`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:457-457`
<!-- vulcan:connections:end -->

## Limitations
Returns `NaN` rather than throwing for a zero or negative propellant budget, and `rpo_group_metric_mean` then silently drops such rows from the average. Values above 100 % are not flagged even though they indicate the vehicle would have run dry.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 456.

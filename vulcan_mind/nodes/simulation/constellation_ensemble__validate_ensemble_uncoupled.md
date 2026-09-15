---
id: simulation.constellation_ensemble__validate_ensemble_uncoupled
label: _validate_ensemble_uncoupled
kind: function
source:
  file: src/simulation/campaigns/constellation_ensemble.jl
  symbol: _validate_ensemble_uncoupled
  lines:
  - 52
  - 52
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: allow_gnc_effectors
  type: Bool
  units: n/a
  required: true
  description: Positional argument `allow_gnc_effectors`.
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
  description: Return value of `_validate_ensemble_uncoupled`.
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

# _validate_ensemble_uncoupled

## Purpose
Guards `run_constellation_ensemble` against configurations whose GNC effectors might couple satellites together. Since each ensemble member is propagated in isolation, any effector that observes or commands another spacecraft (an RPO planner, a coordinated maneuver) would silently produce wrong results, so the check refuses such configurations unless the caller opts in.

## Design & Implementation
If `allow_gnc_effectors` is true the function returns `nothing` immediately. Otherwise it inspects three vectors: `args.guidance_model.guidance_effectors`, `args.navigation_model.navigation_effectors` and `args.control_model.control_effectors`, pushing the field name of each non-empty one into a `coupled_surfaces::Vector{String}`. When that list is empty it returns `nothing`; otherwise it throws an `ArgumentError` whose message joins the offending field names with commas and explains the two remedies: pass `allow_gnc_effectors=true` if every effector acts on a single satellite, or use the monolithic `run_simulation` path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `allow_gnc_effectors` | Bool | n/a | yes | Positional argument `allow_gnc_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_validate_ensemble_uncoupled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · feedback · `src/simulation/campaigns/constellation_ensemble.jl:142-142`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:55-55`
- `callees` → [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
The test is purely presence-based: it cannot distinguish a genuinely single-satellite effector from a coupling one, so any non-empty effector list is rejected by default even when harmless. Conversely, with `allow_gnc_effectors=true` nothing is validated at all. Dynamic effectors in `args.dynamics_model.dynamic_effectors` are not examined even though they too could couple spacecraft.

## Provenance
Mapped from `src/simulation/campaigns/constellation_ensemble.jl` line 52.

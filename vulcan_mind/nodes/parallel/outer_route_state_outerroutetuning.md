---
id: parallel.outer_route_state_outerroutetuning
label: OuterRouteTuning
kind: struct
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: OuterRouteTuning
  lines:
  - 39
  - 39
inputs:
- id: inner_sat_threshold
  type: Int
  units: n/a
  required: false
  description: Field `inner_sat_threshold` (default `8`).
- id: inner_link_threshold
  type: Int
  units: n/a
  required: false
  description: Field `inner_link_threshold` (default `12`).
- id: outer_light_sat_threshold
  type: Int
  units: n/a
  required: false
  description: Field `outer_light_sat_threshold` (default `2`).
- id: outer_light_link_threshold
  type: Int
  units: n/a
  required: false
  description: Field `outer_light_link_threshold` (default `4`).
- id: outer_light_mission_threshold_s
  type: Float64
  units: n/a
  required: false
  description: Field `outer_light_mission_threshold_s` (default `14_400.0`).
- id: spice_constellation_process_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `spice_constellation_process_enabled` (default `true`).
- id: spice_constellation_min_sats
  type: Int
  units: n/a
  required: false
  description: Field `spice_constellation_min_sats` (default `4`).
- id: adaptive_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `adaptive_enabled` (default `true`).
- id: adaptive_min_samples
  type: Int
  units: n/a
  required: false
  description: Field `adaptive_min_samples` (default `2`).
- id: adaptive_exploration_c
  type: Float64
  units: n/a
  required: false
  description: Field `adaptive_exploration_c` (default `1.25`).
- id: failure_penalty_s
  type: Float64
  units: n/a
  required: false
  description: Field `failure_penalty_s` (default `120.0`).
- id: mc_process_min_samples
  type: Int
  units: n/a
  required: false
  description: Field `mc_process_min_samples` (default `16`).
- id: mc_process_min_mission_s
  type: Float64
  units: n/a
  required: false
  description: Field `mc_process_min_mission_s` (default `3600.0`).
- id: process_max_workers
  type: Int
  units: n/a
  required: false
  description: Field `process_max_workers` (default `Sys.CPU_THREADS`).
- id: trace
  type: Bool
  units: n/a
  required: false
  description: Field `trace` (default `false`).
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
  type: OuterRouteTuning
  units: n/a
  description: Constructed `OuterRouteTuning` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# OuterRouteTuning

## Purpose
Immutable bundle of thresholds and adaptive-selection knobs for the outer parallel-route policy. It decides how many satellites or links justify inner threading, when a workload is light enough to run serially, when a SPICE constellation or Monte Carlo batch should go to process workers, and how the adaptive bandit explores.

## Design & Implementation
`Base.@kwdef struct` with defaults: `inner_sat_threshold = 8`, `inner_link_threshold = 12`, `outer_light_sat_threshold = 2`, `outer_light_link_threshold = 4`, `outer_light_mission_threshold_s = 14_400.0` (4 hours), `spice_constellation_min_sats = 4`, `mc_process_min_samples = 16`, `mc_process_min_mission_s = 3600.0`. Adaptive selection uses `adaptive_min_samples = 2` before trusting a route's mean, an upper-confidence exploration constant `adaptive_exploration_c = 1.25`, and `failure_penalty_s = 120.0` seconds added for failed runs. `process_max_workers` defaults to `Sys.CPU_THREADS` because process workers run with `--threads=1` and are limited by physical cores, not `Threads.nthreads()`. `trace::Bool` enables routing diagnostics.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `inner_sat_threshold` | Int | n/a | no | Field `inner_sat_threshold` (default `8`). |
| in | `inner_link_threshold` | Int | n/a | no | Field `inner_link_threshold` (default `12`). |
| in | `outer_light_sat_threshold` | Int | n/a | no | Field `outer_light_sat_threshold` (default `2`). |
| in | `outer_light_link_threshold` | Int | n/a | no | Field `outer_light_link_threshold` (default `4`). |
| in | `outer_light_mission_threshold_s` | Float64 | n/a | no | Field `outer_light_mission_threshold_s` (default `14_400.0`). |
| in | `spice_constellation_process_enabled` | Bool | n/a | no | Field `spice_constellation_process_enabled` (default `true`). |
| in | `spice_constellation_min_sats` | Int | n/a | no | Field `spice_constellation_min_sats` (default `4`). |
| in | `adaptive_enabled` | Bool | n/a | no | Field `adaptive_enabled` (default `true`). |
| in | `adaptive_min_samples` | Int | n/a | no | Field `adaptive_min_samples` (default `2`). |
| in | `adaptive_exploration_c` | Float64 | n/a | no | Field `adaptive_exploration_c` (default `1.25`). |
| in | `failure_penalty_s` | Float64 | n/a | no | Field `failure_penalty_s` (default `120.0`). |
| in | `mc_process_min_samples` | Int | n/a | no | Field `mc_process_min_samples` (default `16`). |
| in | `mc_process_min_mission_s` | Float64 | n/a | no | Field `mc_process_min_mission_s` (default `3600.0`). |
| in | `process_max_workers` | Int | n/a | no | Field `process_max_workers` (default `Sys.CPU_THREADS`). |
| in | `trace` | Bool | n/a | no | Field `trace` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OuterRouteTuning | n/a | — | Constructed `OuterRouteTuning` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.default_outer_route|default_outer_route]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:314-314`
- [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:373-373`
- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:538-538`
- [[parcore.outer_route_metrics_record_outer_route_feedback_bang|record_outer_route_feedback!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_metrics.jl:15-15`
- [[simulation.adaptive_routing__campaign_route_tuning|_campaign_route_tuning]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:119-119`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Defaults are captured at struct construction, so `process_max_workers = Sys.CPU_THREADS` reflects the coordinator's host and is wrong when workers run on a machine with a different core count. No field is range-checked; a negative `adaptive_exploration_c` or a zero `adaptive_min_samples` is accepted silently. Thresholds are absolute counts and seconds, not scaled by per-satellite cost, so a few very expensive satellites can be routed as light.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 39.

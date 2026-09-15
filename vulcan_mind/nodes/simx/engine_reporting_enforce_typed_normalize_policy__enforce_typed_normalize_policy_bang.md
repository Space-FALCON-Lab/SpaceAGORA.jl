---
id: simx.engine_reporting_enforce_typed_normalize_policy__enforce_typed_normalize_policy_bang
label: _enforce_typed_normalize_policy!
kind: function
source:
  file: src/simulation/engine/reporting.jl
  symbol: _enforce_typed_normalize_policy!
  lines:
  - 13
  - 26
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Configuration whose simulation_settings.normalize flag decides whether
    the legacy normalized pipeline was requested.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: policy_decision
  type: Nothing
  units: n/a
  description: Returns nothing when the configuration is acceptable; otherwise throws
    an ArgumentError naming the required change.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# _enforce_typed_normalize_policy!

## Purpose
`_enforce_typed_normalize_policy!` is the gate that keeps the typed propagation path SI-native. It is called at the top of `run_simulation`, before any state is built, and it either lets the run proceed, warns once, or refuses the configuration outright.

## Model & Assumptions
Typed propagation works in metres, seconds and kilograms throughout, so the legacy `SimulationSettings.normalize` flag has no meaning in it: honouring the flag would silently rescale states that downstream effectors interpret as SI. The policy is therefore refusal by default. The escape hatch exists only so legacy transition checks can still be run and is expected to be temporary.

## Design & Implementation
The function returns immediately when `normalize` is false, which is the common case and costs one field read. When the flag is set it consults `_typed_allow_transition_normalize()`, backed by `SPACEAGORA_ALLOW_TYPED_NORMALIZE`. If the transition flag is on it delegates to `_warn_typed_normalize_transition_flag!`, which checks `_typed_normalize_warning_enabled()` and the module-level `_normalize_warning_emitted` reference so the warning is emitted at most once per process rather than once per run, then returns nothing. Otherwise it throws an `ArgumentError` whose message states that typed propagation is SI-native and names both remedies: set `normalize=false`, or set the environment variable for legacy checks only.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `args` | SimulationConfiguration | n/a | yes | Configuration whose simulation_settings.normalize flag decides whether the legacy normalized pipeline was requested. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `policy_decision` | Nothing | n/a | — | Returns nothing when the configuration is acceptable; otherwise throws an ArgumentError naming the required change. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:163-163`

**Downstream**

- `callees` → [[simulation.reporting__warn_typed_normalize_transition_flag_bang|_warn_typed_normalize_transition_flag!]] · `callers` · call · `src/simulation/engine/reporting.jl:18-18`
- `callees` → [[simulation.setup__typed_allow_transition_normalize|_typed_allow_transition_normalize]] · `callers` · call · `src/simulation/engine/reporting.jl:17-17`
<!-- vulcan:connections:end -->

## Limitations
The once-per-process warning latch is a plain `Ref{Bool}` with no lock, so under concurrent campaign workers the warning can be emitted more than once or, in a benign race, be observed as already emitted. The latch is never reset, so a long-lived session that runs many normalized legacy configurations sees the warning only for the first.

## Provenance
Mapped from `src/simulation/engine/reporting.jl:13-26`, with the once-only warning helper at line 1 of the same file and the call site at `src/simulation/engine/execution.jl:163`.

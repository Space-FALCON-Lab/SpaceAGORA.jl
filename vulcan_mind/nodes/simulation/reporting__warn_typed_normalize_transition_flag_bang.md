---
id: simulation.reporting__warn_typed_normalize_transition_flag_bang
label: _warn_typed_normalize_transition_flag!
kind: function
source:
  file: src/simulation/engine/reporting.jl
  symbol: _warn_typed_normalize_transition_flag!
  lines:
  - 1
  - 1
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Nothing
  units: n/a
  description: Return value of `_warn_typed_normalize_transition_flag!`; mutates `args`
    in place. Returns `nothing`.
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

# _warn_typed_normalize_transition_flag!

## Purpose
Emits a single deprecation warning when a typed `run_simulation` is launched with the legacy `SimulationSettings.normalize = true` flag while the transition escape hatch is active, telling the user that typed propagation is SI-native.

## Design & Implementation
Returns `nothing` immediately if `args.simulation_settings.normalize` is false or `_typed_normalize_warning_enabled()` is false. It then checks the module-level `Ref` `_normalize_warning_emitted[]`; if already true it returns, otherwise it sets the flag to `true` before issuing the `@warn`, giving at-most-once semantics per process. The warning text states that propagation is always SI-native in metres, seconds, and kilograms, and advises setting `normalize = false` to silence it. It is invoked from `_enforce_typed_normalize_policy!`, which calls it only on the `_typed_allow_transition_normalize()` branch and otherwise throws `ArgumentError`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_warn_typed_normalize_transition_flag!`; mutates `args` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_reporting_enforce_typed_normalize_policy__enforce_typed_normalize_policy_bang|_enforce_typed_normalize_policy!]] · `callees` → `callers` · call · `src/simulation/engine/reporting.jl:18-18`

**Downstream**

- `callees` → [[simulation.setup__typed_normalize_warning_enabled|_typed_normalize_warning_enabled]] · `callers` · call · `src/simulation/engine/reporting.jl:2-2`
<!-- vulcan:connections:end -->

## Limitations
The at-most-once guard is a process-global `Ref` with no lock, so two threads entering simultaneously can both observe `false` and both warn; conversely the flag is never reset, so a second simulation run in the same session with the same misconfiguration produces no warning at all. The suppression is per-process rather than per-run, which means a long-lived worker process warns only for whichever run happens to be first.

## Provenance
Mapped from `src/simulation/engine/reporting.jl` line 1.

---
id: gnc.propulsive_maneuvers__control_effector_log_enabled
label: _control_effector_log_enabled
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _control_effector_log_enabled
  lines:
  - 10
  - 10
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Bool
  units: n/a
  description: Return value of `_control_effector_log_enabled`.
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

# _control_effector_log_enabled

## Purpose
Decides whether the control-effector layer should emit diagnostic warnings for this simulation run.

## Design & Implementation
Returns `true` outright when `ENV["SPACEAGORA_DEBUG_CONTROL"]` equals the exact string `"1"`. Otherwise it probes the ODE parameter object inside a `try`: `hasproperty(p, :shared_buffers)`, then `hasproperty(p.shared_buffers, :debug_control)`, then `Bool(p.shared_buffers.debug_control[])` dereferencing the `Ref`. Any thrown error is caught and the function returns `false`, so a hand-built parameter object can never crash the logging check.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_control_effector_log_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__control_effector_exception_fallback|_control_effector_exception_fallback]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:327-327`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the literal `"1"` enables the environment override, so `true`, `on` or `yes` are silently ignored, unlike the parsing used elsewhere in the codebase. The bare `catch` hides a genuinely malformed `debug_control` buffer, turning a configuration error into permanently silent diagnostics.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 10.

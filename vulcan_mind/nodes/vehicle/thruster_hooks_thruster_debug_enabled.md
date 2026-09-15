---
id: vehicle.thruster_hooks_thruster_debug_enabled
label: thruster_debug_enabled
kind: function
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: thruster_debug_enabled
  lines:
  - 10
  - 10
inputs:
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
  description: Return value of `thruster_debug_enabled`. Returns `get(ENV, "SPACEAGORA_DEBUG_THRUSTER",
    "0") == "1"`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# thruster_debug_enabled

## Purpose

Reports whether per-firing thruster debug logging is switched on. It reads the environment variable `SPACEAGORA_DEBUG_THRUSTER` and returns `true` only when its value is exactly the string `"1"`, defaulting to `"0"` when the variable is absent.

## Design & Implementation

Declared `@inline` so the check collapses to a cheap string comparison at the call site inside `thrust_calculation_schmitt_trigger!`, where it guards an append-mode `CSV.write` to `thruster_debug.csv`. Keeping the gate as a function rather than a compile-time constant lets a running session toggle logging by mutating `ENV` between simulation runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `thruster_debug_enabled`. Returns `get(ENV, "SPACEAGORA_DEBUG_THRUSTER", "0") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- [[vehicle.thruster_hooks_thrust_calculation_schmitt_trigger_bang|thrust_calculation_schmitt_trigger!]] · `callees` → `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:68-68`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The lookup hits `ENV` on every thruster update rather than caching the result, and only the literal value `"1"` enables logging - `"true"`, `"yes"` or `"ON"` are all treated as disabled. There is no per-vehicle or per-thruster granularity: the switch is global to the process.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl` line 10.

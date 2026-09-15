---
id: simulation.setup__density_without_aero_warning_enabled
label: _density_without_aero_warning_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _density_without_aero_warning_enabled
  lines:
  - 66
  - 66
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
  description: Return value of `_density_without_aero_warning_enabled`. Returns `_engine_env_get("SPACEAGORA_WARN_DENSITY_WITHOUT_AERO",
    "1") == "1"`.
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

# _density_without_aero_warning_enabled

## Purpose
Gate for the diagnostic that warns when a non-vacuum density model is configured without any effector that converts atmosphere into force, so heating-only studies can opt out of the noise.

## Design & Implementation
Returns `_engine_env_get("SPACEAGORA_WARN_DENSITY_WITHOUT_AERO", "1") == "1"`. Enabled by default; the source comment cites the GRAM quickstart as a legitimate density-without-drag case that sets the variable to `0`. Queried once at the top of `_warn_density_without_atmospheric_effector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_density_without_aero_warning_enabled`. Returns `_engine_env_get("SPACEAGORA_WARN_DENSITY_WITHOUT_AERO", "1") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__warn_density_without_atmospheric_effector|_warn_density_without_atmospheric_effector]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:88-88`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:66-66`
<!-- vulcan:connections:end -->

## Limitations
Exact-string comparison means values such as `"true"` disable rather than enable the warning. Suppressing the warning also suppresses it for cases where the omission is a genuine configuration bug.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 66.

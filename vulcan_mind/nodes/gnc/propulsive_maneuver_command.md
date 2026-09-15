---
id: gnc.propulsive_maneuver_command
label: PropulsiveManeuverCommand
kind: function
source:
  file: src/gnc/command_types.jl
  symbol: PropulsiveManeuverCommand
  lines:
  - 5
  - 15
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: CommandTypes namespace exposing the immutable maneuver record to guidance
    and runtime code.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: command
  type: PropulsiveManeuverCommand
  units: n/a
  description: Immutable keyword-constructed maneuver command record shared by guidance,
    control, and runtime state.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
- command
charts:
- gnc
origin: agent
---

# PropulsiveManeuverCommand

## Purpose
`PropulsiveManeuverCommand` is the immutable command record that carries a guidance decision into the runtime state. It separates maneuver validity and source-orbit metadata from the mutable buffers that guidance and control use during propagation.

## Theory & Math
The record contains command metadata rather than a continuous control law. A valid maneuver can be interpreted by downstream thruster guidance as a requested burn or attitude action; an invalid command represents the absence of an actionable maneuver. The physical acceleration is computed later from actuator and mass properties.

## Model & Assumptions
Keyword defaults create an inert command with `valid=false`, allowing a failed or not-yet-computed guidance path to leave an explicit no-command value. Consumers must check validity before treating fields as a burn request. The record is expected to cross guidance, control, and ODE parameter boundaries without mutation.

## Design & Implementation
`command_types.jl` declares and exports the struct near the top of the GNC module. `PropulsiveBurnPlan` extends the command concept with burn timing, thrust, specific impulse, impulse, and propellant fields. Guidance hooks write command values into shared buffers, while runtime types store per-spacecraft command vectors for the RHS and callbacks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | CommandTypes namespace exposing the immutable maneuver record to guidance and runtime code. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `command` | PropulsiveManeuverCommand | n/a | — | Immutable keyword-constructed maneuver command record shared by guidance, control, and runtime state. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_sharedbuffers|SharedBuffers]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:735-735`
- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:122-122`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/command_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The record does not validate burn feasibility, thruster availability, attitude reachability, or propellant margin. A valid flag can be true for a command that later violates actuator limits. Since the record is immutable, correction requires constructing and assigning a replacement rather than updating fields in place.

## Provenance
Mapped from `src/gnc/command_types.jl:5-15`.

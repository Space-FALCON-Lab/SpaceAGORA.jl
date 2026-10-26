# GNC Aerobraking Boundary Contract

## Canonical ownership
- `src/models/gnc/guidance/aerobraking/*` owns the core side of typed energy-depletion guidance (EDG): `typed_edg/` holds the model types, platform services and the forwarders of the two EDG decisions, and `energy_depletion_lifecycle.jl` holds run initialization and passage resets.
- `src/models/gnc/control/` owns the EDG control lifecycle, passage heat accounting and targeting control. Its `aerobraking/targeting_control.jl` applies the companion's control decision to the panels; the numerical calculations stay in the companion.
- The optional `SpaceAGORAEDG.jl` package owns the EDG numerical kernels and is the only EDG implementation. SpaceAGORA no longer carries the older E-EDG and T-EDG strategy routes or their strategy selector.

## Required boundaries
- Guidance files must not include control source files directly.
- Numerical EDG kernels must not be added back to core; `typed_edg/algorithms.jl` forwards only the two decisions, `guidance_decision!` and `control_decision!`. No `src` module imports another EDG kernel name from `EDGAlgorithms` or reaches one through it.
- Propulsive maneuver coupling must cross the typed command boundary (`PropulsiveManeuverCommand`, `AerobrakingControlCommand`) rather than `DynamicEffectors` or direct thruster-state mutation from guidance.

## Canonical interfaces
- `AerobrakingEnergyDepletionConfig`
- `AerobrakingEnergyDepletionState`
- `AerobrakingEnergyDepletionGuidanceModel`
- `AerobrakingEnergyDepletionControlModel`
- `PropulsiveManeuverCommand`
- `AerobrakingControlCommand`
- `SpaceAGORA.EDGServices` (the versioned companion boundary)


## Typed EDG prediction and panel commands

The optional `SpaceAGORAEDG.jl` package owns the typed numerical kernels.
Core's `SimulationModel.EDGAlgorithms`, under
`src/models/gnc/guidance/aerobraking/typed_edg/`, holds two decision forwarders
rather than a second implementation. `SpaceAGORA.EDGServices` is the versioned
boundary. The [EDG package boundary](../src/reference/edg_adapter_boundary.md)
defines the supported calls, units, state ownership, installation and errors.

Guidance and control model types remain in `GuidanceHooks` and `ControlHooks`.
Core owns lifecycle, scheduling, physical integration and panel actuation. Since
contract 2.0.0, and so under SpaceAGORA 0.6's contract 3.0.0, the provider
receives the run-owned EDG state first, then the configuration, copied
observations and call-scoped prediction queries. The
prediction clock is the one the engine reports to that state when the run starts
and at each solve's start (`gnc_clock_started!`); core does not search the
configuration for it.
`ODEParams`, live physical state and mutable vehicle geometry remain inside core.
Shared typed numerical calls forward to the companion.

The guidance decision is also the EDG guidance family of the GNC algorithm
contract (`EDGGuidanceFamily`, in `edg_guidance_family.jl`): `AlgorithmGuidance`
can run it, and trusted core code validates each decision before installing it
in the EDG state.

The control hook requests `control_decision!` at its existing control times and
applies the returned `EDGControlDecision.command` only when `apply` is true.
Prediction uses `config.controlled_panel_links`, which may differ from the
control effector's links. An invalid index produces no command application or
environment query. A `:solved` status records execution, not scientific acceptance.

Numerical decisions may update the run's EDG plan, cache and diagnostic state;
they do not move panels or physical heat integrals. Exceptions still propagate
without rollback: a decision can update mode before an atmosphere error, and
completed telemetry can precede an actuator error. Control-first scheduling,
passage bookkeeping and panel application order remain core responsibilities.

Targeting retains its trajectory-based, wind-relative atmosphere query. Heat-load
prediction retains its altitude-only, zero-latitude/longitude query. These are
different scientific assumptions, not interchangeable implementations. Packaging
does not change query order, callback timing, constraints or acceptance criteria.

The GNC boundary and source-completeness gates reject numerical implementations
added back to core. Focused ownership tests cover command/application parity,
invalid-index laziness, switch endpoints, panel-selection mismatch and query
conventions. Package installation and ownership checks do not replace integrated
mission or scientific validation.

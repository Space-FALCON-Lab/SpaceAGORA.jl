# GNC aerobraking boundary contract

Energy-depletion guidance (EDG) is split between core and the optional
`SpaceAGORAEDG.jl` package. The package owns the numerical kernels and is the
only EDG implementation; core owns the run lifecycle, scheduling, physical
integration and panel actuation. The older E-EDG and T-EDG strategy routes and
their strategy selector have been removed from core.

## Ownership

- `src/models/gnc/guidance/aerobraking/` holds the core side of typed EDG.
  `typed_edg/` contains the model types, the platform services and the
  forwarders of the two EDG decisions; `energy_depletion_lifecycle.jl` handles
  run initialization and passage resets; `edg_guidance_family.jl` defines the
  EDG guidance family of the GNC algorithm contract.
- `src/models/gnc/control/aerobraking/` holds the EDG control lifecycle
  (`energy_depletion_lifecycle.jl`), passage heat accounting
  (`edg_heat_accounting.jl`) and targeting control (`targeting_control.jl`).
  `targeting_control.jl` applies the companion's control decision to the
  panels; the numerical work stays in the companion.
- `SimulationModel.EDGAlgorithms`, in `typed_edg/algorithms.jl`, holds two
  forwarders, `guidance_decision!` and `control_decision!`, rather than a second
  implementation.
- `SpaceAGORA.EDGServices` (`src/models/gnc/companions/edg/services.jl`) is the
  versioned boundary the companion builds against. The
  [EDG package boundary](../src/reference/edg_adapter_boundary.md) defines its
  supported calls, units, state ownership, installation and errors.

The canonical interfaces are `AerobrakingEnergyDepletionConfig`,
`AerobrakingEnergyDepletionState`, `AerobrakingEnergyDepletionGuidanceModel`,
`AerobrakingEnergyDepletionControlModel`, `PropulsiveManeuverCommand`,
`AerobrakingControlCommand` and `SpaceAGORA.EDGServices`. The guidance and
control model types are defined inside the `GuidanceHooks` and `ControlHooks`
modules, which include the files above.

## Required boundaries

- Guidance files do not include control source files.
- Numerical EDG kernels are not added back to core. `typed_edg/algorithms.jl`
  forwards only `guidance_decision!` and `control_decision!`, and no `src`
  module imports any other EDG kernel name from `EDGAlgorithms` or reaches one
  through it.
- Propulsive maneuver coupling crosses the typed command boundary
  (`PropulsiveManeuverCommand`, `AerobrakingControlCommand`). Guidance does not
  go through `DynamicEffectors` or change thruster state directly.

`test/gates/ci_gnc_aerobraking_boundary_gate.jl` and
`test/gates/ci_src_completeness_contract_gate.jl` refuse numerical code added
back to core, using the checks in `test/contracts/edg_ownership_checks.jl`.

## Decisions and prediction

Under the current service contract (3.0.0, `CONTRACT_VERSION` in
`typed_edg/adapter_queries.jl`) the provider receives the run's EDG state first,
then the configuration, a copied observation and a prediction context that is
valid only for that call. The prediction clock is the one the engine reports to
that state when the run starts and again at the start of each solve
(`gnc_clock_started!`); core does not search the configuration for it.
`ODEParams`, the live physical state and the mutable vehicle geometry stay
inside core: the companion sees copied values through the prediction context.

The guidance decision is also the EDG family of the GNC algorithm contract
(`EDGGuidanceFamily`). `AlgorithmGuidance` can run it, and trusted core code
validates each decision (`validate_edg_guidance_decision`) before installing it
in the EDG state.

The control hook calls `control_decision!` at its existing control times and
applies the returned `EDGControlDecision.command` only when its `apply` field is
true. Some predictions inside the companion use `config.controlled_panel_links`,
which can differ from the control effector's links. An out-of-range spacecraft
index returns a decision with `apply` false and makes no environment query. A
`:solved` switch action records that the switch solve ran, not that its result
was scientifically accepted.

Numerical decisions may update the run's EDG plan, cache and diagnostic state.
They do not move panels or change physical heat integrals. Exceptions propagate
without rollback: a decision can update the guidance mode before an atmosphere
error, and completed telemetry can precede an actuator error. Control-first
scheduling, passage bookkeeping and the order of panel application stay core
responsibilities.

Targeting keeps its trajectory-based, wind-relative atmosphere query
(`_edg_targeting_prediction_environment`). Heat-load prediction keeps its
altitude-only query at zero latitude and longitude
(`_edg_sample_prediction_atmosphere`). These are different scientific
assumptions, not interchangeable implementations. Packaging changes neither the
query order nor callback timing, constraints or acceptance criteria.

## Tests

`test/gnc/aerobraking/energy_depletion_ownership_tests.jl` covers command and
panel-application parity, lazy handling of an invalid index, switch endpoints,
panel-selection mismatch and the query conventions above. Installation and
ownership checks do not replace an integrated mission or scientific validation.

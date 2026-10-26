# GNC and planning ownership

Guidance decides what a vehicle should do. Planning produces a reference that
can be checked and tracked. Control turns that reference into actuator
commands, and navigation provides the hook for a sensor or estimator update.
These parts share state through explicit contracts. Their source folders under
`src/models/gnc/` are a layout, not Julia modules: a file loads into the module
that includes it, whichever folder it sits in.

## Planning and its adapters

`src/models/gnc/planning/` holds the planning sources in three groups. `rpo/`
has the RPO planner contract (`contract/`), target-relative geometry
(`geometry/`), the direct baseline planner, the RRT kernels, sampling,
retiming, path metrics, the planner comparison (`comparison/`) and the
planner-neutral lifecycle. `robot_arm/` has the quintic joint references, the
joint-reference family and the arm path geometry. `common/` has the path and
tree operations those applications share.

Core's side of the companions is in `src/models/gnc/companions/`. `hypr/` holds
the core-owned configuration, availability checks, the `HYPRServices` facade
and the `HYPRRPOPlanner` adapter that the optional HYPR package plugs into, and
`edg/` holds `EDGServices`, the facade of the optional SpaceAGORAEDG package.

HYPR's search implementation stays in its own repository. A `HYPRServices` name
in core identifies an extension point or a service that core supplies to that
implementation; it is not a second copy of the search. Module and type
identities do not change when files move between folders, so the bindings in
`GuidanceHooks`, `RobotArmPlanning`, `HYPRUtils` and `NavigationHooks` keep the
names existing code uses. The [HYPR service contract](../src/reference/hypr_services.md)
governs the external boundary.

RPO has one integration path. The public planner interface connects a planner
to `RPOPlannerLifecycle`, which validates each returned reference before
installing it in the buffer it shares with one controller. HYPR plans through
it as `HYPRRPOPlanner`; the HYPR-specific `RPOGuidanceModel`, which read its
HYPR configuration directly, was removed in 0.5.0. See the
[RPO pilot guide](../src/user/rpo_planner_pilot.md).

The `navigation/` directory holds the navigation extension hook. The RPO
reference and collision geometry is in `planning/rpo/geometry/` but loads into
`NavigationHooks`, where its bindings have always been. That geometry describes
the planner's world and answers clearance queries; it is not a navigation
estimator.

## Comparison tracking and reporting

Planner selection and batch orchestration stay with the planning sources.
Comparison tracking is implemented in the private
`src/interfaces/analysis/planning/rpo_comparison_tracking.jl`, which loads as
`SimulationModel.RPOComparisonTracking` after guidance and control. It imports
the controller operations from `ControlHooks`, the actuator type from
`ThrusterModels`, and the reference, geometry and plant operations it uses.
The HCW and two-body tracking loops keep their separate behavior.

The six tracking generic functions are declared in `GuidanceHooks`, and
`RPOComparisonTracking` adds their only methods. Qualified calls, the
`SimulationModel` bindings and serialized function references therefore still
resolve to the same defining names. Tracking settings, planner configuration
and the two-body plant type are defined where they were before. This internal
module adds nothing to the root public API.

EDG pairing validation finds the EDG control model by dispatch.
`GuidanceHooks._is_edg_control` returns `false` by default, and `ControlHooks`
adds the method for `AerobrakingEnergyDepletionControlModel`. Guidance never
looks up `ControlHooks`, or any other module, through its parent module at run
time; it gets the EDG decisions from `EDGAlgorithms`, which loads before
`GuidanceHooks` and is imported by it. `test/unit/gnc/planning_ownership_tests.jl`
fails on such a lookup.

Reporting on completed comparisons is in
`src/interfaces/analysis/planning/rpo_comparison_reports.jl`. Its internal
`RPOComparisonReports` module writes the result CSVs and prepares the metrics
and figure data. It imports from `GuidanceHooks` only the report functions'
original identities and the planner label helpers; it does not choose a
planner or change a reference, a controller or an acceptance result.

`planning/rpo/comparison/reporting_compat.jl` declares those report functions
in `GuidanceHooks` and then loads the reporting module, which adds their only
methods. Existing `GuidanceHooks` and `SimulationModel` imports, qualified
calls and serialized function references keep resolving to the same names.
The module is an implementation boundary, not public API.

The Plotly figures come from the PlotlyJS package extension
(`ext/SpaceAGORAPlotlyJSExt.jl`), through
`src/interfaces/analysis/planning/rpo_comparison_plots.jl`. Without PlotlyJS
loaded, a figure function, or a comparison writer asked for figures, throws
`PlotlyJSUnavailableError`; a writer called with `write_plotly_outputs=false`
and `write_failed_path_outputs=false` writes the CSV without figures.

## EDG

The independent, experimental `SpaceAGORAEDG.jl` companion owns the EDG
numerical kernels. Core keeps the public model types, the run-owned state and
its resets between passages, the physical integration, the queries, the
outputs and the panel actuation. `guidance/aerobraking/typed_edg/` holds these
core models and services. Its `algorithms.jl` forwards the two decisions,
`guidance_decision!` and `control_decision!`, through the private
`SimulationModel.EDGServices` module in `typed_edg/services.jl`, which copies
the observation and lends a prediction context for the call, to the
`SpaceAGORAEDGExt` extension; the extension accepts only the names listed in
the public `SpaceAGORA.EDGServices.IMPLEMENTED_FUNCTIONS`. Control calls the
companion's `control_decision!` and applies the command to the panels; it
holds no EDG calculation and binds no EDG helper.

SpaceAGORA no longer carries the older E-EDG and T-EDG strategy routes: its
port of the original energy depletion guidance (Falcone and Putnam, JGCD
45(4), 2022) and of the earlier targeting route. The typed path and the
companion implement the generalized formulation (Simha, Yu and Falcone,
SciTech 2026) and make no claim to reproduce the 2022 paper exactly. To run the
old routes, check out tag `parallel-paper-freeze-2026-03-05`, the last tag with
the routes wired into the simulation path, or tag `v1.0.0`, whose older example
scripts configure them.

The [EDG service contract](../src/reference/edg_adapter_boundary.md) defines
the supported boundary, and the [EDG guide](../src/user/edg_guidance.md)
covers installation and the limits of the experimental scope. Passing the
extraction and interface conformance checks is not a broader scientific
qualification.

## Robotics

Robot hardware stays with the vehicle models (`src/models/vehicle/robotics/`),
and the compliant multibody equations with the dynamics
(`src/models/dynamics/multibody/`). MuJoCo works through the engine's
[external-propagation contract](../src/reference/external_propagation_services.md),
which governs who propagates the state, not planning. The `packages/robotics/`
folder that holds the HYPR and SpaceAGORAMuJoCo checkouts is only a grouping;
it moves none of these responsibilities to either companion.

## Companion boundaries

[Companion boundaries](companion_boundaries.md) lists, for HYPR, SpaceAGORAEDG
and SpaceAGORAMuJoCo, the private SpaceAGORA paths their tests and research
scripts still use, the rules for moving them to supported calls, and the checks
to run when you change a boundary.

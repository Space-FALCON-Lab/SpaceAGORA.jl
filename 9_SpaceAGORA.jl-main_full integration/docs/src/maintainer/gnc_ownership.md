# GNC and planning ownership

Guidance decides what a vehicle should do. Planning produces a reference that
can be checked and tracked. Control turns that reference into actuator commands,
while navigation provides the hook for a sensor or estimator update. These
responsibilities share state through explicit contracts; their source folders do
not create new Julia modules.

## Planning and its adapters

`src/models/gnc/planning/` collects the planning sources in four groups. `rpo/`
contains target-relative geometry, the direct baseline, RRT kernels, sampling,
retiming, metrics and the planner-neutral lifecycle. `robot_arm/` contains quintic
joint references and arm path geometry. `common/` contains path and tree
operations shared by those applications. `adapters/hypr/` contains the core-owned
configuration and compatibility surface used by the optional HYPR implementation.

HYPR's search implementation remains in its companion repository. A
`HYPRServices` name in core may identify an extension point or a service supplied
to that implementation; it does not imply another copy of the search algorithm.
Existing module and type identities are retained, including historical
`GuidanceHooks`, `RobotArmPlanning`, `HYPRUtils` and `NavigationHooks` bindings.
The [HYPR contract](hypr_services.md) governs the external boundary.

RPO has two integration paths. The public planner interface connects a planner
to `RPOPlannerLifecycle`, which validates returned references before installing
them in the buffer shared with one controller. The existing `RPOGuidanceModel`
uses its HYPR configuration directly and remains available for existing studies.
Their configuration and lifecycle assumptions must not be interchanged merely
because both produce an `RPOPlan`. See the [public planner guide](../user/rpo_planner_pilot.md).

The `navigation/` directory now contains the navigation extension hook. RPO
reference and collision geometry lives in `planning/rpo/geometry/`; its defining
`NavigationHooks` bindings remain unchanged. That geometry describes the planner's
world and clearance queries. It is not, by itself, a navigation estimator.

## Reporting is an analysis responsibility

Planner comparison execution and tracking stay with the planning sources.
Reporting completed results belongs to
`src/interfaces/analysis/planning/rpo_comparison_reports.jl`. Its internal
`RPOComparisonReports` module owns CSV generation, metric presentation and Plotly
output. It imports only the historical report function identities and the planner
label helpers from `GuidanceHooks`; it does not choose a planner or change a
reference, controller or acceptance result.

The `reporting_compat.jl` file declares the old function objects before loading
the analysis module. The analysis module adds their sole implementations, so
existing `GuidanceHooks` and `SimulationModel` imports, qualified helper access,
and serialized function references continue to resolve to the same defining
names. The namespace is an implementation boundary, not a new public API. This
change does not make plotting an optional package dependency.

## EDG and legacy bridges

EDG remains part of core GNC. Typed EDG prediction and decision calculations are
in `guidance/aerobraking/typed_edg/`. Its guidance model selects modes and targets;
its control model applies panel commands. Their shared state, atmosphere-passage
resets and refusal rules remain unchanged. The short files in `control/` that
import EDG calculations are compatibility bindings, not additional algorithms.

The older `e_edg/` and `t_edg/` strategies and `control/aerobraking/` execution
paths remain separate. `internal/bridge_helpers.jl` supplies their legacy runtime
context conversions. This organization does not retire either path or certify
their scientific equivalence. The [EDG boundary plan](edg_adapter_boundary.md)
describes the further validation needed before a separate adapter can be adopted.

Robot hardware remains under vehicle models, and compliant multibody equations
remain under dynamics. MuJoCo uses the engine's
[external-propagation contract](external_propagation_services.md), which governs
propagation ownership rather than planning. A Robotics grouping for companion
checkouts does not transfer these core responsibilities to either companion.

# Guidance, navigation, control and planning

This is a source-owner map, not a new Julia namespace hierarchy. Existing module,
type and function identities remain stable, including qualified compatibility names.

Guidance and control keep their models in the same application subfolders:
`aerobraking/`, `landing/`, `maneuvers/` and `rpo/`, and control also
`robot_arm/`, `attitude/` and `interlinks/`.

| Directory | Responsibility |
| --- | --- |
| `contracts/` | The general GNC contract: command types, the algorithm contract and its simulation binding, the runtime declarations and the lifecycle hooks. These are domain contracts, not the top-level CLI/file-interface group. |
| `guidance/` | Simulation guidance models (`models.jl`) and update hooks (`hooks.jl`): typed EDG in `aerobraking/`, Apollo descent in `landing/`, thruster guidance in `maneuvers/`, and in `rpo/` the hooks that consume planned references. |
| `control/` | Tracking controllers, held actuator commands and simulation effectors (`hooks.jl`): EDG heat accounting and targeting control in `aerobraking/` (the EDG calculations are the companion's; control calls its decision), Apollo descent in `landing/`, propulsive maneuvers and the burn-plan family in `maneuvers/`, LQ-MPC in `rpo/`, the arm controller in `robot_arm/`, the magnetic momentum manager in `attitude/`, interlink scheduling in `interlinks/`. |
| `navigation/` | Sensor/estimator extension hook (`hooks.jl`). Historical `NavigationHooks` geometry identities load their implementations from `planning/rpo/geometry/`. |
| `planning/rpo/` | The RPO planner contract (`contract/`), direct baseline, RPO geometry, reusable sampling/retiming/metrics, RRT and comparison planner kernels, and the neutral planner lifecycle. |
| `planning/robot_arm/` | Quintic joint-space references and reusable arm path geometry. Vehicle construction and physical dynamics have separate owners. |
| `planning/common/` | Path geometry and RRT tree operations shared by RPO and robot-arm planning. |
| `companions/` | Core-owned adapters of the optional companion packages. `edg/`: the `EDGServices` facade of SpaceAGORAEDG. `hypr/`: HYPR configuration, availability and versioned services, compatibility entry points and simulation integration; HYPR search execution remains in the external HYPR repository. |

Completed planner-comparison CSVs and plots are implemented in
`src/interfaces/analysis/planning/rpo_comparison_reports.jl`, inside the internal
`RPOComparisonReports` module. GNC's `reporting_compat.jl` retains the historical
`GuidanceHooks` function objects; analysis imports those identities and supplies
their methods. No second report implementation or forwarding method is installed.

The planner-neutral `RPOPlannerLifecycle` validates results before updating the
shared controller reference; it is the only RPO guidance path (the HYPR-specific
`RPOGuidanceModel` was removed in 0.5.0). For typed EDG,
`guidance/aerobraking/typed_edg/` owns core model types, services and compatibility
dispatch; the optional SpaceAGORAEDG package owns the numerical kernels. Core
retains lifecycle and panel actuation, and its shared typed calculations forward to
the companion. The older E-EDG/T-EDG strategy routes were removed; the companion is
the only EDG implementation. See the [EDG package boundary](../../../docs/src/reference/edg_adapter_boundary.md)
for ownership and supported calls.

Interlink scoring and terminal matching live in `control/interlinks/interlink_scheduling.jl`.
They choose actuation from a candidate
graph; the graph and terminal metadata belong to vehicle models. Accepted-step
callback lifecycle belongs to simulation orchestration, and the central laser
force belongs to coupled dynamics. These responsibilities use the existing
source domains.

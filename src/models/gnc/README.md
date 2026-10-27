# Guidance, navigation, control and planning

This folder groups the GNC code by role. The folders are not Julia modules:
each file's module is set by its `module` line and by the file that includes
it, and the older qualified names (for example `GuidanceHooks.<name>`) still
work.

Guidance and control use the same application subfolders: `aerobraking/`,
`landing/` and `maneuvers/`, plus `rpo/` in both, and `robot_arm/`,
`attitude/` and `interlinks/` in control only.

| Folder | What it holds |
| --- | --- |
| `contracts/` | The general GNC contract: command types (`CommandTypes`), the algorithm contract (`GNCAlgorithms`) and its simulation binding (`AlgorithmLifecycle`), the run-time declarations (`GNCRuntime`) and the lifecycle hooks (`SimulationLifecycle`). |
| `guidance/` | The guidance models (`models.jl`, module `GuidanceModels`) and update hooks (`hooks.jl`, `GuidanceHooks`): typed EDG in `aerobraking/`, Apollo descent in `landing/`, thruster guidance in `maneuvers/`, and in `rpo/` the plan buffer (`RPOPlan`, `RPOPlanBuffer`) that the RPO planner lifecycle writes and the RPO controller reads. |
| `control/` | Tracking controllers, held actuator commands and the control hooks (`hooks.jl`, `ControlHooks`): EDG heat accounting and targeting control in `aerobraking/` (the EDG calculations come from the companion; control calls its decision), Apollo descent in `landing/`, propulsive maneuvers and the burn-plan family in `maneuvers/`, LQ-MPC in `rpo/`, the arm controller in `robot_arm/`, the magnetic momentum manager in `attitude/`, and interlink scheduling in `interlinks/`. |
| `navigation/` | The sensor and estimator hook (`hooks.jl`, `NavigationHooks`). The older RPO geometry names in `NavigationHooks` load their code from `planning/rpo/geometry/`. |
| `planning/rpo/` | The RPO planner contract (`contract/`), the direct baseline planner, RPO geometry, sampling, retiming and metrics, RRT and the comparison planners, and the planner lifecycle (`RPOPlannerLifecycle`). |
| `planning/robot_arm/` | Quintic joint-space references and arm path geometry, and the joint-reference family. Building the arm and its dynamics are done elsewhere. |
| `planning/common/` | Path geometry and RRT tree operations shared by RPO and robot-arm planning. |
| `companions/` | SpaceAGORA's side of the optional companions. `edg/services.jl` is the public `SpaceAGORA.EDGServices` facade for SpaceAGORAEDG. `hypr/` holds HYPR's configuration, availability checks, the `HYPRServices` facade, entry points and simulation integration; the HYPR search itself runs in the HYPR repository. |

Three algorithm families are loaded by `src/SpaceAGORA.jl` itself rather than
by a hooks module: the burn-plan family (`control/maneuvers/burn_planning.jl`,
`BurnPlanning`), the joint-reference family
(`planning/robot_arm/joint_reference_family.jl`, `JointReferencePlanning`) and
the EDG guidance family (`guidance/aerobraking/edg_guidance_family.jl`,
`EDGGuidanceFamilies`).

## Planner-comparison reports

`planning/rpo/comparison/reporting_compat.jl` declares the report functions
inside `GuidanceHooks`, so their old qualified names keep working, and then
includes `src/interfaces/analysis/planning/rpo_comparison_reports.jl`, which
defines the `RPOComparisonReports` module there and supplies the methods. That
file includes `rpo_comparison_plots.jl` for the figures, which need the PlotlyJS
extension (`ext/SpaceAGORAPlotlyJSExt.jl`). `rpo_comparison_tracking.jl`, in the
same folder, supplies the comparison's tracking methods. There is one
implementation of each report; nothing forwards to a second copy.

## RPO and EDG

`RPOPlannerLifecycle` checks each planner result before it updates the shared
controller reference, and it is the only RPO guidance path; the HYPR-specific
`RPOGuidanceModel` was removed in 0.5.0.

For typed EDG, `guidance/aerobraking/typed_edg/` holds SpaceAGORA's model
types, the internal `SimulationModel.EDGServices` queries (atmosphere, frames
and state projections) and the dispatch that keeps older calls working. The
optional SpaceAGORAEDG package makes the two decisions, through the public
`SpaceAGORA.EDGServices` facade in `companions/edg/`. SpaceAGORA keeps the run
lifecycle and the panel actuation. The older E-EDG/T-EDG strategy routes were
removed, so the companion is the only EDG implementation. See the
[EDG package boundary](../../../docs/src/reference/edg_adapter_boundary.md) for
ownership and supported calls.

## Interlinks

Interlink scoring and terminal matching are in
`control/interlinks/interlink_scheduling.jl` (`SchedulingPolicies`). They choose
which links to use from a candidate graph; the graph and the terminal data
belong to the vehicle models. The callback that runs the schedule after each
accepted step belongs to simulation orchestration, and the laser force on the
spacecraft belongs to the dynamics code.

# Guidance, navigation, control and planning

This is a source-owner map, not a new Julia namespace hierarchy. Existing module,
type and function identities remain stable, including qualified compatibility names.

| Directory | Responsibility |
| --- | --- |
| `guidance/` | Simulation guidance models and update hooks; aerobraking/EDG, landing and propulsive decisions. RPO hooks consume planned references. |
| `control/` | Tracking controllers, held actuator commands and simulation effectors. The small typed EDG calculation files retain compatibility bindings. |
| `navigation/` | Sensor/estimator extension hook. Historical `NavigationHooks` geometry identities load their implementations from `planning/rpo/geometry/`. |
| `planning/rpo/` | Direct baseline, RPO geometry, reusable sampling/retiming/metrics, RRT and comparison planner kernels, and the neutral planner lifecycle. |
| `planning/robot_arm/` | Quintic joint-space references and reusable arm path geometry. Vehicle construction and physical dynamics have separate owners. |
| `planning/common/` | Path geometry and RRT tree operations shared by RPO and robot-arm planning. |
| `planning/adapters/hypr/` | Core-owned HYPR configuration, availability and versioned services, compatibility entry points and simulation integration. HYPR search execution remains in the external HYPR repository. |
| `interfaces/` | Planner request/result/validation contracts and simulation lifecycle extension points. These are domain contracts, not the top-level CLI/file-interface group. |
| `internal/` | Legacy aerobraking state/context bridges shared by guidance and control. |

Completed planner-comparison CSVs and plots are implemented in
`src/interfaces/analysis/planning/rpo_comparison_reports.jl`, inside the internal
`RPOComparisonReports` module. GNC's `reporting_compat.jl` retains the historical
`GuidanceHooks` function objects; analysis imports those identities and supplies
their methods. No second report implementation or forwarding method is installed.

The planner-neutral `RPOPlannerLifecycle` validates results before updating the
shared controller reference. The existing `RPOGuidanceModel` remains a separate
HYPR-specific entry path used by existing studies. Typed EDG remains in core:
`guidance/aerobraking/typed_edg/` owns its calculations, while the older E-EDG/T-EDG
paths and their control bridges remain available. Folder changes do not establish
numerical equivalence between those scientific formulations.

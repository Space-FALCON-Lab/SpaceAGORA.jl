# Vulcan Run Contract

**Run ID:** `run-1782181174229-yx08xn` (also in env var `VULCAN_RUN_ID`)

The Vulcan MCP server is available at `http://127.0.0.1:7711` (also in env var `VULCAN_MCP_URL`).

You are running as a non-interactive agent. Use these Vulcan MCP tools:

- `start_run(run_id)` — call once when you begin (optional but recommended)
- `list_sheets()`, `get_sheet(sheet_id)`, `query_nodes(...)` — read the model
- `read_component(component_id)`, `update_component(...)` — read/write per-component docs
- `update_task_column(task_id, column, reason)` — move Kanban cards as you work
- `write_design_notes(agent_name, content, mode)` — persistent memory across runs
- `report_progress(run_id, message)` — short status updates (visible to user)

**REQUIRED before you exit:** call `complete_run` exactly once:

```json
complete_run({
  run_id: "run-1782181174229-yx08xn",
  status: "success" | "partial" | "failed",
  summary: "<short markdown describing what you changed>",
  modified_sheet_ids: [...],     // optional
  open_questions: [...]          // optional; each becomes a Kanban clarification card
})
```

If you do not call `complete_run`, Vulcan treats this run as partial/failed and shows the user a generic error.

---

---

# Task: Implement: ads_sensor_models.jl

StarTracker/Gyro/Magnetometer/SunSensor measure() + noise.

## Mission Context
SpaceAGORA.Jl
Mission: Align SpaceAGORA.Jl with the repository implementation and keep the model current.
Compile intent: Align the implementation with the current system model and close repository-to-model gaps.

## Source Change
add-node: New module block added to sheet "ADS — Software Architecture"

## Requirements You Must Satisfy
No requirements are linked to this task. You have freedom of implementation approach.

## V&V Tests You Must Pass
No V&V tests are linked. Run the standard test suite; no specific test targets required.

## Operational Context
No operational flow sheets linked.

## Physical Interface Constraints
No physical interface constraints. Implementation approach is unconstrained by hardware.

## Allowed Components and Files
### CommandTypes (src/gnc)
Files: src/gnc/command_types.jl, src/gnc/control/aerobraking/constraint_tracking.jl, src/gnc/control/aerobraking/control_commands.jl, src/gnc/control/aerobraking/tracking_executor.jl, src/gnc/control/attitude_control_models.jl, src/gnc/control/control_hooks.jl, src/gnc/control/propulsive_maneuvers.jl, src/gnc/guidance/aerobraking/common/closed_form_solution.jl, src/gnc/guidance/aerobraking/common/heat_rate_models.jl, src/gnc/guidance/aerobraking/dispatcher.jl, src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl, src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl, src/gnc/guidance/aerobraking/e_edg/second_switch_solver.jl, src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl, src/gnc/guidance/aerobraking/interfaces.jl, src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl, src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl, src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl, src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl, src/gnc/guidance/guidance_hooks.jl, src/gnc/guidance/guidance_models.jl, src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl, src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl, src/gnc/internal/bridge_helpers.jl, src/gnc/navigation/navigation_hooks.jl
Interfaces: PropulsiveManeuverCommand(export/out), PropulsiveBurnPlan(export/out), AerobrakingControlCommand(export/out), calcControlForceTorque(export/out), calcControlEffect!(export/out), calcControlMassFlowRate(export/out), calcGuidanceEffect!(export/out), AbstractAerobrakingStrategy(export/out), EEdgStrategy(export/out), TEdgStrategy(export/out), AerobrakingGuidanceInput(export/out), AerobrakingGuidanceOutput(export/out), compute_aerobraking_guidance(export/out), dispatch_aerobraking_guidance(export/out), AerobrakingCampaignPropulsiveManeuverGuidanceModel(export/out), calcNavigationEffect!(export/out)

### NavigationHooks (src/gnc/navigation)
Files: src/gnc/navigation/navigation_hooks.jl, src/gnc/navigation/ads_sensor_models.jl, src/gnc/navigation/ads_mekf.jl, src/gnc/navigation/ads_mode_manager.jl, src/gnc/navigation/ads_navigation_effector.jl
Interfaces: calcNavigationEffect!(export/out)

## Constraint Manifest
Your constraint manifest is at: C:/Users/Robbie/falcon/SpaceAGORA.jl/.vulcan/projects/spaceagora-jl/runs/run-1782181174229-yx08xn/manifest.json
It contains: allowedFiles, frozenSignatures, linkedRequirements, requiredVerifiers.
Read it before making any changes.

## Behavioral Rules (ENFORCED by MCP tools, not advisory)
1. Use propose_file_patch — not direct file edits. Vulcan validates and applies patches.
2. Use propose_component_patch — to update ComponentDoc metadata.
3. If a file is outside allowedFiles: call request_scope_expansion. Do NOT edit it directly.
4. If completing the task requires a new interface, port, or component: call declare_infeasible.
5. Call run_verifier for ALL required verifiers before complete_run.
6. Use write_run_note continuously: reasoning, math, issues, decisions, test results.

## Escalation
Scope elevation (need access to another file/component):
  request_scope_expansion({ requested_files, requested_components, rationale, blocking_constraint, proposed_resolution, risks })

Infeasibility (model must change for this to be implementable):
  declare_infeasible({ reason, required_changes: [{kind, description, affected_model_elements}], proposed_resolution, risks })

## Required Verifiers
Run these before complete_run: interface-diff, model-consistency
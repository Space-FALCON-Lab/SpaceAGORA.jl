# Missions

Mission definitions, operational sequencing and schedules. Simulation
configuration, the solver and campaigns are in `../simulation/`.

`mission_workflow.jl` defines `MissionPlan`, named observation phases
(`MissionPhase`), mission events (`MissionEvent`) and `run_mission`.
`mission_segments.jl` defines operational segments, which move from one to the
next on a deadline, after a duration or when a condition holds, all within one
continuous solve.

The three behave differently at run time. Observation phases only label the
results; they switch no models and add no solver stops. Segments run their
declared entry and exit actions through solver callbacks that
`mission_segments.jl` builds. Each mission event builds its callback with the
factory you give it. The [mission workflow guide](../../../docs/src/user/mission_workflows.md)
covers the supported workflow and its limits.

`operations/maneuver_plans.jl` holds reusable Odyssey firing schedules, with
their recorded values and direction conventions, and the conversion to typed
maneuver commands. The package does not load it: `examples/AGORA_Odyssey.jl`
and the integration test harness include it explicitly.

GNC applies the actuator commands; the optional SpaceAGORAEDG package computes
the typed EDG decisions.

# Mission responsibilities

This directory owns mission definitions, operational sequencing and schedules.
Simulation configuration, solver execution and campaigns remain in `../simulation/`.

`mission_workflow.jl` owns `MissionPlan`, named observation phases, event factories
and `run_mission`. `mission_segments.jl` owns operational segments and their
deadline, duration or condition transitions within one continuous solve.
Observation phases label results; segments execute declared entry and exit
actions. See the [mission workflow guide](../../../docs/src/user/mission_workflows.md)
for the supported workflow and limits.

`operations/maneuver_plans.jl` contains reusable Odyssey firing schedules and the
adapter to typed maneuver commands. `examples/AGORA_Odyssey.jl` includes it
explicitly; it is an active scenario dependency, even though the package entry
point does not load it. The schedules retain their recorded values and direction
conventions.

GNC owns guidance and control integration and applies actuator commands; the
optional SpaceAGORAEDG package owns typed EDG numerical prediction and decisions.

Mission configuration describes the intended run. Runtime phase transitions and
events are executed by the simulation callbacks; a folder reorganization must
not change their triggers, order or solver interaction.

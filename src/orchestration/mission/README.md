# Mission responsibilities

This directory owns mission-specific schedules and policy selection. Simulation
configuration, solver execution and campaigns remain in `../simulation/`.

`operations/maneuver_plans.jl` contains reusable Odyssey firing schedules and the
adapter to typed maneuver commands. `examples/AGORA_Odyssey.jl` includes it
explicitly; it is an active scenario dependency, even though the package entry
point does not load it. The schedules retain their recorded values and direction
conventions.

`operations/aerobraking_policy/` owns the default strategy selector and per-pass
schedule. It selects the strategy requested by the mission; EDG numerical
prediction and the actuator commands remain in GNC guidance and control.
`DRLPolicyAdapterStub` is an explicit unimplemented extension point. It throws
when selected and must not be described as a trained or validated policy.

Mission configuration describes the intended run. Runtime phase transitions and
events are executed by the simulation callbacks; a folder reorganization must
not change their triggers, order or solver interaction.

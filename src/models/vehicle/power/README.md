# Power

The spacecraft power budget.

| File | Holds |
|---|---|
| `power_components.jl` | the records, owned by `Components`: `SolarArray`, `Battery` (with its run state: full, empty, shedding tier) and `PowerLoad` |
| `power_system.jl` | `PowerSystems`: `add_power!`, generation, the bus load, the actuators' power by shedding tier, the battery's rate, the shedding mask and the run-state updates |

The engine integrates the battery's stored energy (`battery_wh`, a state only
when a battery is mounted) in a pass after the right-hand side
(`src/orchestration/simulation/engine/dynamics_rhs.jl`) and applies the
shedding mask to the mounted actuators' levels. The solver modes that split the
right-hand side (`:gravity_backbone_split`, `:split_imex`, `:multirate`) refuse
a run with a battery. The battery's state switches at events and step ends in
`src/orchestration/simulation/callbacks/power_callbacks.jl`, which also writes
the power columns.

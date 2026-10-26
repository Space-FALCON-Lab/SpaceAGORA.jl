# Power

The spacecraft power budget.

| File | Holds |
|---|---|
| `power_components.jl` | the records (owned by `Components`): `SolarArray`, `Battery` (with its run state: full, empty, shedding tier), `PowerLoad` |
| `power_system.jl` | `PowerSystems`: `add_power!`, generation, the bus load, the actuators' power by shedding tier, the battery's rate, the shedding mask and the run-state updates |

The engine integrates the battery's stored energy (`battery_wh`, in the state
only when a battery is mounted) in a pass after every right-hand-side route
(`src/orchestration/simulation/engine/dynamics_rhs.jl`), applies the shedding
mask to the mounted actuators' levels, and switches the battery's state at
events and step ends (`src/orchestration/simulation/callbacks/power_callbacks.jl`),
which also writes the power columns.

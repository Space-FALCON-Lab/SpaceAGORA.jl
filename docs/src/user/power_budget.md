# Powering a Spacecraft

A spacecraft can carry a power budget: solar arrays that generate, a battery
that stores, and loads that draw. The mounted actuators draw the power they
report through `actuator_consumption` (laser terminals, and reaction wheels
and torque rods given power figures). When a battery is mounted, its stored
energy is part of the integrated state, and when it runs empty the engine
sheds actuator load until the battery has recharged.

## The parts

| Part | Mounted on | What it does |
|---|---|---|
| `SolarArray(; normal, area_m2, efficiency, max_power_w=Inf)` | any link | generates `area_m2 * efficiency * S * max(0, cos θ) * ν` W, where `S` is the solar flux, `θ` the angle between `normal` (link frame) and the Sun, and `ν` the lit fraction of the Sun's disk |
| `Battery(; capacity_wh, initial_soc, charge_efficiency, discharge_efficiency, max_charge_w, max_discharge_w, recovery_soc=0.05)` | the root link, at most one | stores the surplus and covers a shortfall |
| `PowerLoad(; power_w, name)` | any link | a constant load, never shed |

Mount them with `add_power!(spacecraft, link, part)`. An array on a hinged
panel turns with the panel. Without `initial_soc`, the battery starts at the
spacecraft's `battery_energy_index` (a percentage).

## What happens during a run

- **Generation.** The Sun's direction comes from the run's ephemerides: SPICE,
  or, with `SimpleEphemeridesModel`, a low-precision model for Earth (other
  bodies need SPICE). The solar flux uses the central body's mean distance
  from the Sun (1 AU for Earth and the Moon); the eclipse model is the one of
  solar radiation pressure. A `PowerReceiver` adds the power its incoming power
  beams deliver (see [Interlinks](interlinks.md#Power-links)); a receiver alone
  counts as a power part.
- **Charging.** A surplus charges the battery at `charge_efficiency`, up to
  `max_charge_w`; a full battery curtails the rest.
- **Discharging.** A shortfall discharges it over `discharge_efficiency`, up
  to `max_discharge_w`.
- **Shedding.** When the battery is empty, the engine sheds actuator load in
  tiers until what is left fits what the arrays generate: laser terminals
  (and the radios that transmit on data links) first, then torque rods, then the reaction wheels' motors (a shed wheel puts
  no torque on the vehicle). Shedding ends once the battery has recharged to
  `recovery_soc` of its capacity. Constant loads are never shed; while the
  battery is empty and they exceed the generation, the net power stays
  negative.

The battery's full and empty states switch at events, when its stored energy
reaches its capacity or zero; the shedding tier is updated at the end of each
integration step. The split solver modes (`:split_imex`, `:multirate`,
`:gravity_backbone_split`) refuse a spacecraft with a battery. The opt-in
`:multirate_mri` route integrates battery and interlink power in its fast
component and supports the same boundary events. A checkpointed
run keeps these states from one segment to the next, and a run resumed from a
checkpoint file takes them from the file, which holds the stored energy and
whether the battery is empty and recharging or full.

## Add power parts to a spacecraft

This configuration fragment builds the vehicle and mounts its power parts.
It does not yet define an orbit or run a simulation.

```julia
using SpaceAGORA, StaticArrays
bus = Link(root=true, m=10.0, dims=MVector(0.1, 0.2, 0.3))
spacecraft = SpacecraftModel(links=[bus], root=bus, inertia_tensor=bus.inertia)
add_power!(spacecraft, bus, SolarArray(normal=SVector(1.0, 0.0, 0.0), area_m2=0.06, efficiency=0.3))
add_power!(spacecraft, bus, Battery(capacity_wh=5.0, initial_soc=0.5))
add_power!(spacecraft, bus, PowerLoad(power_w=5.0, name="avionics"))
add_actuator!(spacecraft, bus, LaserTerminal(boresight=SVector(1.0, 0.0, 0.0), max_power_w=2.0, efficiency=0.2))
```

For a complete runnable scenario, use
[`Earth_CubeSat_Power_Budget.jl`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/0487e8226107592b19221618f0f21f500cf2aacf/examples/Earth_CubeSat_Power_Budget.jl)
with the environment described in the [Examples Catalog](examples_catalog.md).
It flies this vehicle for two orbits
with its laser on: the battery fills in sunlight, empties in eclipse, the
laser is shed until the Sun returns, and the script checks each of these.

The results carry the [power columns](output_columns.md#Power-budget):
stored energy, state of charge, generation, load, shed power and net power.

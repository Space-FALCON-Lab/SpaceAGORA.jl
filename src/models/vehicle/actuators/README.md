# Actuators

Hardware mounted on a spacecraft link that acts on the vehicle when commanded.
Each kind describes its hardware once; the engine derives the force, torque,
propellant use, power draw and result columns from the actuators mounted on
the vehicle, and samples the environment the mounted kinds ask for.
The extensibility page and `templates/actuator_template.jl` show how to add a
kind.

| File | Holds |
|---|---|
| `actuators.jl` | the contract (`Actuators`): the methods every kind implements, the evaluation context and the environment a kind asks for (`ActuatorRequirements`), the controller hooks `actuation_mode` and `actuator_levels`, and `allocate_actuators` |
| `mounting.jl` | mounting and the mounted set (`ActuatorMounting`): `add_actuator!`, the mounting order, `allocation_matrix`, `mounted_wrench`, the applied levels and power behind the result columns |
| `thruster.jl` | the thruster kind: the `Thruster` record (its `direction` is the thrust direction; `firing_dynamics` opts a mounted thruster into the firing dynamics), its contract methods and `six_axis_thrusters` |
| `laser.jl` | the laser-terminal kind: the `LaserTerminal` record (level: the transmitted optical power; the vehicle recoils with `level / c`), its contract methods, and the link parameters it gives the interlink model |
| `fault.jl` | the `ActuatorFault` record (failed, stuck, degraded or biased between two times) and the spacecraft's `FaultSet`, owned by `Components` |
| `faults.jl` | `ActuatorFaults`: `add_fault!`, the faulted levels the mounted path uses, and the fault switching the engine runs at the faults' window ends |
| `hinge.jl` | the hinge kind: the `PanelHinge` record (level: the commanded angle of the link it turns; angle and rate limits), applied at control ticks |
| `cmg.jl` | the control-moment-gyro kind: the `ControlMomentGyro` record (level: the gimbal rate; rotor momentum, gimbal and spin axes, rate limit, gimbal stops), `cmg_momentum`, `cmg_torque_axis`, `cmg_angle_saturate`; the `CMGInnerGimbal` record of a double-gimbal CMG and `cmg_torque_axes`; mounted CMGs and inner gimbals live in the root link's wheel assembly and their gimbal angles in the `cmg_angles` state (the outer ones first) |
| `gimbal.jl` | the thruster-gimbal kind: the `ThrusterGimbal` record (level: the gimbal angle of the thruster it turns, about an axis in the link frame; angle and rate limits), applied at control ticks; the thruster's `gimbal_rotation` carries the turn |
| `thruster_models.jl` | `BaseThrusterModel` (a burn schedule per spacecraft; with `actuation=:mounted` it fires a mounted engine and points it with the mounted wheels, the control code being in `gnc/control/maneuvers/propulsive_maneuvers.jl`) and `SixAxisThrusterModel` (six fixed jets) |
| `thruster_models_module.jl` | the `ThrusterModels` owner module |
| `thruster_hooks.jl` | thruster firing dynamics over a link's jets: Schmitt trigger, impulse integration, ramp factor; the engine applies them to mounted thrusters with `firing_dynamics` at control ticks |
| `thruster_allocation.jl` | the RPO six-axis force allocation and its wrench |
| `reaction_wheel.jl` | the reaction-wheel kind: the `ReactionWheel` record (level: the torque on the body about its spin axis; torque and momentum limits, unlimited by default; base and motor power with `wheel_power`), its contract methods, and the root link's `ReactionWheelAssembly` whose momenta the engine integrates |
| `magnetorquer.jl` | `Magnet`, a permanent dipole the torque-rod effector applies, and the torque-rod kind: the `Magnetorquer` record (level: the commanded dipole along its axis; resistive power with the square of the dipole), its contract methods |

Records stay owned by `Components` (`src/models/vehicle/spacecraft/components.jl`),
which includes the kind files, so serialized spacecraft keep their type
identity.

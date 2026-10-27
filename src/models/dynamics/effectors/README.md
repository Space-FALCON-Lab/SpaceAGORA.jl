# Forces and torques on a spacecraft

The code here turns physical models into the forces and torques the dynamics
engine integrates, given a spacecraft's state and the environment sampled for
it. Sampling the atmosphere, ephemerides and fields is environment work; this
folder applies those samples to the spacecraft's mass, attitude, geometry and
hardware.

`force_torque_models.jl` defines the `DynamicEffectors` module. It includes one
file per family from `force_torque_models/`, each its own submodule:

| File | Submodule |
|---|---|
| `gravity_effectors.jl` | `GravityEffectors`: point-mass, J2 and gravity-gradient models |
| `aerodynamic_effectors.jl` | `AerodynamicEffectors`: the aerodynamic coefficient models (from `aerodynamic_wrench_models.jl`) and the mesh aerodynamic surrogate (from `aerodynamic_mesh_surrogate.jl`) |
| `perturbation_effectors.jl` | `PerturbationEffectors`: harmonics, third bodies, radiation, magnetics and the LVLH attitude cascade |
| `guidance_models.jl` | `GuidanceModels`: re-exports the propulsive-maneuver guidance models |
| `robot_arm_reaction_effector.jl` | `RobotArmReactionEffectors`: arm joint loads fed back into the spacecraft |
| `reaction_wheel_momentum.jl` | `ReactionWheelMomentum`: wheels with prescribed speed histories |
| `laser_force_effectors.jl` | `LaserForceEffectors`: deprecated laser-link force names that forward to each link's `link_force` |
| `plume_surface_interaction.jl` | `PlumeSurfaceInteraction`: a phenomenological descent-plume model (surface pressure footprint, erosion diagnostic, near-ground thrust augmentation) |

`force_torque_models/gravitational.jl`, included last, says which effectors
model gravity (`is_gravitational`). The accelerometer callback uses it to
remove gravity from the specific force, and the engine uses it to check the
gravity in each spacecraft's effector scope.

`perturbations.jl` is an ordered list of includes from `perturbations/`, all
evaluated in `PerturbationEffectors`: gravity (harmonics, third bodies),
radiation (pressure, eclipse), magnetics and `control/lvlh_cascade.jl`. The
LVLH cascade is a continuous body-torque law evaluated with the equations of
motion, not a sampled GNC callback.

The files in `gravity/` define the spacecraft side of the analytic gravity
models (model types, run-time force helpers, wrenches, the gravity-gradient
effector) inside `DynamicEffectors.GravityEffectors`. Their field kernels are
in `../../environment/gravity/`, and `../../environment/gravity/gravity_models.jl`
includes both sets of files in a fixed order.

Keep one compiled force kernel per model, called by both the scalar and the
batch routes, so that the two routes give the same results.

Thread eligibility and the choice of multibody workers belong to
`../../../orchestration/parallel/`. The older aerodynamic helper names in
`AerodynamicEffectors` are aliases of that shared policy, so existing callers
reach the same functions and the same per-solve cache.

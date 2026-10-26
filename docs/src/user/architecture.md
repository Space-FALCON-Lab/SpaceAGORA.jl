# How a simulation fits together

This page explains how a configuration becomes a simulation, how to choose a
spacecraft model and which routes are supported for adding behavior. For a
runnable example, start with the [Quickstart](quickstart.md).

## One picture of a run

A run starts when a `SimulationConfiguration` is passed to `run_simulation`.
Setup checks the configuration and builds the initial state and the models.
During the solve the engine evaluates forces and torques, advances the state
and runs the scheduled and event callbacks. Control forces and torques enter
the right-hand side of the equations of motion, while navigation, guidance and
control-command updates happen in callbacks. Values are saved as the solve
proceeds and written to the output files.

![Run flow: scenario configuration contains the spacecraft collection; individual spacecraft contain rigid, articulated or compliant assemblies. Environment, spacecraft and GNC supply the repeated dynamics and solver loop; callbacks update commands and save values for outputs.](../assets/architecture_run_flow.svg)

[Download the vector architecture figure (PDF).](../assets/architecture_run_flow.pdf)

The arrows show how data flows at run time. [The Integrated State](integrated_state.md)
describes the quantities the solver advances, and
[Simulation Outputs](outputs.md) describes the saved values and files.

## Choose the parts of a run

| Part | Your choice | Guide |
| --- | --- | --- |
| Spacecraft | Bodies, links, joints, mass properties, actuators and initial conditions | [Simulation Configuration](simulation_configuration.md) |
| Environment | Planet, atmosphere, gravity, ephemerides and other selected loads | [Atmosphere Models](atmosphere_models.md), [GRAMSuite Setup](gramsuite_setup.md) |
| Guidance and control | Configured effectors and supported optional planners | [RPO Planner Pilot](rpo_planner_pilot.md), [HYPR Planning](hypr_planning.md), [EDG Guidance](edg_guidance.md) |
| Propagation | Time span, tolerances, integration method and stopping conditions | [Solver Configuration](solver_configuration.md), [Stopping on a Condition](stop_conditions.md) |
| Execution and outputs | One run or a campaign, process workers, sampling and output files | [Parallel Execution](parallel_execution.md), [Simulation Outputs](outputs.md) |

The no-GRAM route is enough for a first run. GRAM, SPICE and the optional
packages need their own preparation, and only when the models you choose use
them; see [Installation and Environment](installation_environment.md).

## Spacecraft assembly and constellation composition

A spacecraft is built from components, and a scenario holds spacecraft. A
`SpacecraftModel` holds its links (with the actuators mounted on them), its
joints, the root body, the mass properties and the initial condition. The
scenario's `DynamicsModel` holds the collection of spacecraft and the dynamic
effectors that act on them. That collection can be a single vehicle, a
servicing pair or a constellation. On its own it does not assign tasks or
provide a communication network.

### Rigid and articulated spacecraft

A spacecraft becomes articulated when one of its `Joint`s has `joint_type`
`:hinge`, `:slide` or `:ball`; `:fixed` joints merge their links into the
parent body. The joint coordinates are integrated beside the root state as
`joint_q` and `joint_qd`, and the root is a free-floating base. By default the
existing non-gravity loads act on the root body, computed from the configured
geometry, and gravity is evaluated for each body at its own center of mass.
With `SimulationSettings(articulated_live_pose_loads=true)`, aerodynamics and
facet SRP are evaluated for each link at its live pose, and the gravity
gradient becomes a per-body torque. A spacecraft with no joint other than
`:fixed` does not use the articulated path at all. See
[Joint types](articulated_compliant.md#Joint-types-and-articulated-spacecraft).

### Compliant assemblies

A compliant multibody model represents flexibility as discrete bodies joined
by springs and dampers; you specify the bodies' masses and inertias and the
joints' stiffness and damping. The supported cases are the articulated and compliant assemblies
documented in [Articulated and Compliant Vehicles](articulated_compliant.md);
not every flexible surface or contact interaction is coupled and validated. A
mesh describes surface geometry, while a multibody model describes motion and
loads.

A `CompliantAttachment` mounts a `CompliantMultibodyModel` on a spacecraft's
root, on a link merged into it, or on a moving articulated link. Its bodies are
integrated beside the spacecraft and exchange forces and torques with that
link in both directions. The compliant robot-arm coupling is a separate
supported path. In that path each arm link feels the run's position-only
gravity at its own inertial position, so the joint springs carry the loads that
differ from the bus's free fall. See
[Compliant attachments](articulated_compliant.md#Compliant-attachments) and
the [Compliant Multibody API](../generated/compliant_multibody_api.md). The
[API migration note](../public_api_policy.md#Compliant-multibody-and-robot-arm-names)
explains the compatibility names; they do not add MuJoCo contact coupling.

### Several spacecraft

An arm attached to a station belongs to the station's assembly, while a
free-flying inspector is a spacecraft of its own. Assigning tasks to them and
coordinating them is mission-level work, above the construction of any single
vehicle.

When spacecraft interact, simulate them together in one run.
`run_constellation_ensemble` instead runs each spacecraft as an independent
member, so it refuses an interlink model, and by default it also refuses
configuration-level guidance, navigation and control effectors, which could
couple spacecraft. Its opt-in, `allow_gnc_effectors=true`, is for effectors
that each act on every spacecraft on its own (`spacecraft_scope`
`AllSpacecraft()`); it does not add communication between vehicles. GNC
declared on a `SpacecraftModel` travels with that spacecraft and needs no
opt-in. An independent ensemble and a coordinated constellation run therefore
make different promises.

## Add a model or use an optional package

For a new force, torque, atmosphere or control model, start from the
documented root API and [Extensibility](../extensibility.md); the
[force and torque walkthrough](custom_effector.md) shows a complete example.
Through the supported hooks your model takes part in the existing run while
keeping the units, frames and solver behavior each hook describes.

Larger optional integrations go through versioned service modules:

- HYPR provides the optional planner behind `HYPRRPOPlanner`. SpaceAGORA keeps
  the planner contract, the direct baseline planner and the shared path
  geometry, and supplies configuration, shared kinematics and dynamics through
  the [HYPR services](../reference/hypr_services.md).
- External integrators, including the MuJoCo companion, use the
  [External Propagation Services](../reference/external_propagation_services.md),
  which state which spacecraft each integrator owns and which couplings are
  supported.
- The experimental SpaceAGORAEDG package owns the guidance kernels; SpaceAGORA
  keeps the physical dynamics, the callbacks and the run lifecycle, and connects
  to it through the [EDG services](../reference/edg_adapter_boundary.md).
- The optional visualizer reads saved scene data and does not change how the
  simulation propagates; see [Interactive Visualization](visualization.md).

A source checkout of an optional package does not make it a required
dependency of SpaceAGORA. Use the
[installation commands](installation_environment.md#Optional-package-commands)
for the integration you need. A successful installation or interface check
does not mean the results are scientifically valid for a particular
application.

The [Public API and compatibility](../public_api_policy.md) reference tells you
whether a name is supported. Internal modules can be imported, but they are not
an extension contract. A checkout also contains `experimental/`: typed
scaffolds that `using SpaceAGORA` does not load, with no stability guarantee.

The repository layout, the contributor checklist and the full map of source
owners are kept in the repository, in `CONTRIBUTING.md` and the contributor
guide it links to, not on this site.

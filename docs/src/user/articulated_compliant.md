# Articulated and Compliant Vehicles

Use this page when a spacecraft has moving parts: links connected by joints,
such as a robot arm or deployable panels, or attachments connected by springs
and dampers.

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [The Integrated State](integrated_state.md)
- [Compliant Multibody API](../generated/compliant_multibody_api.md)

## Joint types and articulated spacecraft

A spacecraft is **articulated** if and only if at least one of its `Joint`s is not
`:fixed`. Articulated spacecraft integrate their joint coordinates beside the bus state;
every other spacecraft, including one whose joints are all `:fixed`, runs exactly as before
(no tree is built and nothing is validated). Run `examples/Articulated_Panels_Demo.jl` for a
working case.

### Declaring joints

```julia
bus = Link(root=true, m=620.0, dims=MVector(2.05, 2.05, 2.8))
panel = Link(m=10.0, dims=MVector(0.01, 2.85, 1.0), r=MVector(0.0, 2.45, 0.0))   # COM in the bus frame
add_joint!(sc, Joint(bus, SVector(0.0, 1.025, 0.0), panel, SVector(0.0, -1.425, 0.0);
    joint_type=:hinge, axis=[1.0, 0.0, 0.0], stiffness=60.0, damping=4.0, initial_q=0.2))
```

`Joint(link1, p1ᵇ, link2, p2ᵇ; joint_type=...)` connects the parent `link1` to the child `link2`
at a joint point given in the parent frame (`p1ᵇ`) and in the child frame (`p2ᵇ`). Joint
coordinate 0 is the configured geometry: the links' `r` (COM, bus frame) and `q` (attitude
relative to the bus, scalar-last). For non-fixed joints `p1ᵇ` and `p2ᵇ` must map to the same
point in that geometry (1e-9 m); fixed joints merge at the configured geometry as given.

| `joint_type` | Coordinates | `axis` | `stiffness`, `damping` | `rest` | `initial_q`, `initial_qd` |
|---|---|---|---|---|---|
| `:fixed` (default) | none; the child is merged into the parent body | unused | unused | unused | unused |
| `:hinge` | 1 angle (rad) about `axis` | required, parent frame, normalized | scalars (N m/rad, N m s/rad) | angle (default 0) | scalars (default 0) |
| `:slide` | 1 displacement (m) along `axis` | required, parent frame, normalized | scalars (N/m, N s/m) | displacement (default 0) | scalars (default 0) |
| `:ball` | scalar-last quaternion (3 rate components) | unused | scalar or 3x3 PSD matrix | quaternion (default identity) | quaternion (default identity), parent-frame angular velocity 3-vector (default 0) |

Gains must be `>= 0`. Springs and dampers act in joint space: `τ = -k (q - rest) - c q̇`; for a
ball joint the spring uses the axis-angle of `rest⁻¹ ⊗ q` and the damper the relative angular
velocity (the 3x3 stiffness is exact for an isotropic `k I`, first-order for an anisotropic one).
Existing joints keep working and are `:fixed`; the legacy `Kx/Kt/Cx/Ct` fields are unused by this
feature. Links whose joint is `:fixed` move rigidly with their parent body.

### Conventions and state

- The root body is the root link plus every link merged into it through `:fixed` joints plus
  the propellant. The engine's `pos` and `vel` are the **root composite center of mass**; the
  system center of mass is saved as `sc{i}_system_com_*`.
- Quaternions are scalar-last and are the body-to-inertial rotation, like the bus `q`.
- `joint_q` and `joint_qd` are appended to the spacecraft's state, listing the non-fixed
  joints in `spacecraft.joints` order. A ball joint has four coordinates and three rates, and
  its rate is the angular velocity of the child relative to its parent in the parent link
  (`Joint.link1`) frame. Ball coordinates, rates, rest and joint-space gains keep this frame
  even when fixed joints merge the parent link into a body whose head has a different attitude;
  output `joint_q` and `joint_qd` are not converted to the composite body frame.
- Saved outputs: `joint_q`, `joint_qd`, `articulated_link_pose` (every link's inertial COM and
  attitude from the joint kinematics) and `system_com`; see [Simulation Outputs](outputs.md).
- Checkpoint and resume work: the joint state is part of the checkpointed state.

### Mass and inertia

- `mass` in the state keeps its meaning: total spacecraft mass, dry plus propellant, so mass flow
  and `sc1_mass` are unchanged. The root body mass at run time is the state mass minus the
  moving bodies' masses (which never change).
- The link masses and inertias are used. The spacecraft's `inertia_tensor` is **ignored** in
  articulated mode, and `dry_mass` must equal the sum of the link masses (checked at setup).
- Propellant is a point mass at the root composite COM: it adds mass and no inertia and never
  moves the COM, so the root inertia stays the configured composite.

### Loads in v1 (default)

Every non-gravity load (aerodynamics, SRP, thrusters, control torques, third-body gravity) is
evaluated by the existing effectors from the rigid configured geometry and applied to the root
body. Gravity from the position-only gravity models (point mass, J2, harmonics) is evaluated per
body at its own COM, so the inter-body part of the gravity-gradient effect of the layout, including
joint motion, emerges from the dynamics. The default mode cannot apply a gravity-gradient torque
itself, so a gravity-gradient request on an articulated spacecraft (`gravity_gradient=true` on a
gravity effector, or a `GravityGradientTorqueModel`) is **refused at setup** with an `ArgumentError`
that names the effector. Set `SimulationSettings(articulated_live_pose_loads=true)` (below) to apply per-body
gravity gradient from the link inertias, or remove the request. Rigid spacecraft are unaffected. This is what
runs unless the switch below is turned on.

### Live-pose loads (opt-in)

```julia
SimulationSettings(articulated_live_pose_loads=true)    # default false
```

With the switch on, aerodynamics and facet SRP act on each link where the link actually is, and
the gravity gradient becomes per body. Everything is evaluated from one joint-kinematics pass per
RHS call; the `Link` objects (`r`, `q`) are never written.

- **Aerodynamics** (`AerodynamicCoefficientfM`, `AerodynamicCoefficientConstant`,
  `AerodynamicCoefficientNoBallisticFlight`): each link gets its own flow angles from its live
  attitude and its own airspeed from its live COM velocity, which includes the `omega x r` term of
  the body that carries it. That term is where aerodynamic damping of a rotating or swinging
  appendage comes from. The force acts at the link COM; the link's `cop_offset_b` gives the torque
  about the link COM. Density, temperature and wind are the spacecraft-level sample unless the
  model's `per_link_atmosphere` option is set, which queries non-root links at their own position
  as before. The saved drag, lift and cross vectors are the sums over the links.
- **Facet SRP** (`FacetSolarRadiationPressureModel`): facet normals and centers of pressure use the
  live link attitude, and the Sun distance and the eclipse fraction are evaluated per link at the
  link COM, so one panel can be in shadow while another is lit.
- **Where the loads go:** each link's force and torque are applied to the dynamic body that carries
  the link (a link merged by a `:fixed` joint rides with its parent body), as a force through the
  link COM plus the torque about the body COM, `tau_link + (R_body c_link) x F_link`. This also
  removes the v1 error of taking the bus-origin lever about the root composite COM.
- **Gravity gradient:** any request, either `gravity_gradient=true` on `ConstantGravityModel`,
  `InverseSquaredGravityModel` or `InverseSquaredJ2GravityModel`, or a `GravityGradientTorqueModel`,
  becomes the per-body torque `3 mu/r^3 r x (R_b I_b R_b^T) r` from the tree's body inertia `I_b`
  about each body COM. The configured `inertia_tensor` version is not applied as well, so nothing is
  counted twice, and a gravity effector's flag is no longer dropped. Together with the per-body
  point-mass gravity this gives the full gradient torque of the whole structure. The per-body torque
  exists only with the switch; without it the same request is refused (see above). No request, no torque.
- **Everything else** (thrusters, control torques, cannonball SRP, third-body gravity, magnets,
  eddy damping, the mesh aerodynamic surrogate, heating) keeps the base application from v1.
  Compliant attachments stay gravity-only.

Limits and refusals with the switch on:

| Case | Behavior |
|---|---|
| Panel-angle control of a link on a moving body | refused at setup, as without the switch |
| Panel-angle control of a `:fixed` link of the root body | works; the live link attitude follows the commanded `link.q` |
| Whole-vehicle (`root_only`) `AerodynamicCoefficientMeshSurrogate` | refused at setup: the fit bakes in the configured geometry |
| Rigid spacecraft | never enter this path |
| Link mass, COM and inertia | frozen at the configured geometry (unchanged) |

For a locked configuration (very stiff joints) the live-pose run reproduces the rigid run except for
the `omega x r` airspeed term, of relative size `|omega| |r| / |v|`, and a hinged rod in a circular
orbit librates at `sqrt(3) n`, which needs the per-body gradient.

### Solver advice and supported routes

| Route | Support |
|---|---|
| Single-RHS first-order solver modes `:tsit5`, `:auto_stiff`, `:rodas5p`, `:dp8` | supported |
| `:split_imex`, `:multirate`, `:multirate_mri`, `:symplectic`, `:gravity_backbone_split` | refused with an `ArgumentError` |
| Serial, `satellite_batch` and per-satellite RHS routes | supported (one workspace per spacecraft) |
| Flat constellation queue | automatically rerouted to the per-satellite route; forcing `SPACEAGORA_RHS_EXECUTION_MODE=flat` is refused |
| Process routes and `isolate_state` copies | supported: runtime data is rebuilt from the configuration inside each run |
| `orientation_sim=false`, reaction wheels, a robot-arm effector on the same spacecraft | refused |
| Kinematic panel-angle control (`SolarPanelAngleOfAttackControlModel`) | works on links of the root body exactly as before; an error on links of a moving body until joint motors exist |
| Mixing articulated and rigid spacecraft in one run | refused (the state holds equal-sized blocks); articulated spacecraft must share link and joint counts |

Joint dynamics can be stiff: for a stiff hinge keep `dt_max_orbit` below a fraction of the hinge
period, and use tight tolerances; `:dp8` was used for the conservation checks. The joint
coordinates use the quaternion and angular-rate tolerances.

!!! note
    Motors and joint limits are not part of this release. Per-link loads are opt-in; see "Live-pose loads" above.

## Compliant attachments

A `CompliantAttachment` mounts a `CompliantMultibodyModel` (a panel mesh or another
flexible appendage represented by discrete bodies) on a spacecraft link. The attachment bodies move under their own
compliant joint springs, dampers and actuators and exchange forces and torques with the
link in both directions. `examples/Cloth_Panel_Attachment_Demo.jl` runs the four-panel compliant
deployment of `Solar_Panel_Cloth_Deployment_Demo.jl` inside `run_simulation` this way.

```julia
build = build_rectangular_compliant_grid(3, 4; anchor_index=1)       # or any CompliantTopologyBuild / model
att = CompliantAttachment(;
    model=build,                       # a CompliantMultibodyModel, or a build (its state is the initial state)
    link=bus,                          # root, a fixed-merged link, or a link of a moving body
    mount_point=(0.5, 0.0, 0.0),       # mount frame origin in the link frame (m)
    mount_quaternion=(0.0, 0.0, 0.0, 1.0),   # mount frame orientation relative to the link frame
    joint_actuators=actuators,         # CompliantJointActuator list
    rest_schedule=(out, t) -> (out[1] = rest_quaternion(t); nothing),   # optional
)
sc = SpacecraftModel(; links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=ic, attachments=[att])
```

### Frames and conventions

- The **mount frame** sits at `mount_point` in the link frame with orientation
  `mount_quaternion` relative to it. The link frame origin is the link center of mass (the
  spacecraft position for the root link; the link's `r` and `q` relative to the bus otherwise)
  and its axes are the link body axes. In an articulated spacecraft the link frame follows the
  link's dynamic body.
- The model's joints with `parent == 0` attach to the mount frame instead of a fixed base: the
  model's `base_position` and `base_quaternion` are **ignored**. At least one such joint is
  required. A topology build's own state, given for its base pose, is converted to mount-frame
  coordinates; an explicit `initial_state` is already in the mount frame (13 numbers per body:
  position, scalar-last quaternion, velocity, body-frame angular velocity, relative to a mount at
  rest). With neither, the model starts at its rest geometry (a spanning tree from the mount joints;
  pass a build's state for meshes with closed loops).
- The run state (`att_r`, `att_q`, `att_v`, `att_ω`) holds each body's position and velocity RELATIVE to the
  spacecraft's `pos` and `vel` (inertial axes) and its absolute attitude and body rate; the initial state is the
  mount-frame state moved by the mount's pose and rigid-body velocity about the bus. The saved `attachment_pose`
  and `system_com` are inertial (`pos + att_r`). (The robot-arm state `arm_*` still holds absolute positions: a known follow-up.)
- `rest_schedule(out::Vector{SVector{4,Float64}}, t)` writes every joint's rest quaternion for time `t`
  (seconds since the start) into `out`; `nothing` keeps each joint's own rest quaternion. It is called
  on every RHS evaluation with the stage time, so a smoothstep deployment is followed exactly (a
  standalone stepper that freezes the rest per step lags by half a step), and the run is
  allocation-free when the schedule is.

### How the loads reach the spacecraft

The joint loads and the gravity difference are computed from relative quantities. The bus (or articulated base)
acceleration, which includes the reactions, is then known, and each body's relative acceleration is
`F_i/m_i + g(r_base + r_rel) - g(r_base) + (g(r_base) - a_base)`: gravity is evaluated at each body's absolute position
but only its difference from the base gravity enters, so the large common term never reaches the relative acceleration.

- The compliant joint springs, dampers and actuators use the same math as `compliant_joint_loads`
  (`compliant_joint_loads_in_place!` is the allocation-free variant the engine calls).
- The reaction on the mount, a force at the mount point plus a torque, goes to the link it is mounted
  on. **Rigid spacecraft:** the force is added to the bus force (inertial) and the torque, including
  the `r_mount x F` lever, to the bus torque (body frame) before the bus equations. **Articulated
  spacecraft:** the body kinematics are computed first, then the attachment loads, and the reaction enters
  `articulated_dynamics!` as a per-body wrench (`body_force_world`, `body_torque_world`, about the body
  COM) that is mapped to the generalized forces through each body's Jacobian columns, so the hinge,
  slide and ball coordinates and the root respond to the attachment.
- **Loads in v1:** gravity acts on every attachment body at its own position (the position-only
  gravity effectors). No aerodynamics, SRP, thermal or other effector acts on attachment bodies.

### Mass bookkeeping

Attachment bodies are not links: they are not in `links`, `dry_mass` or the state `mass`, and carry
their own masses (`attachment_total_mass(sc)`). The system total is the spacecraft plus the attachments;
the saved `system_com` includes them. A body modeled as a `:fixed` link instead adds to `dry_mass`.
Thrust, mass flow and every other effector see the spacecraft mass only.

### Solver advice and supported routes

| Route | Support |
|---|---|
| `:tsit5`, `:auto_stiff`, `:rodas5p`, `:dp8`, rigid or articulated spacecraft | supported |
| `:split_imex`, `:multirate`, `:multirate_mri`, `:symplectic`, `:gravity_backbone_split` | refused with an `ArgumentError` |
| Serial and per-satellite RHS routes, `isolate_state` copies, checkpoint and resume | supported |
| Flat constellation queue | automatically rerouted to the per-satellite route; forcing `SPACEAGORA_RHS_EXECUTION_MODE=flat` is refused |
| `orientation_sim=false`, a robot-arm effector on the same spacecraft | refused |
| Attachment link not in the spacecraft, model without a joint to the mount | refused |
| Mixing spacecraft with different attachment body counts (or with and without attachments) in one run | refused (equal-sized state blocks) |

Compliant multibody meshes can be stiff. Use `:auto_stiff` or `:rodas5p` when the fastest attachment
mode (about `sqrt(k/m)` for the translational springs and `sqrt(k_rot/I)` for the rotational ones,
plus the damping rate `c/m`) is much faster than the orbit and attitude dynamics of interest; the
explicit modes are fine for mild meshes such as the demo (about 30 rad/s). The attitude and
angular-rate tolerances apply to `att_q` and `att_ω`. Because the attachment state is relative to the bus, spring forces carry no roundoff from orbital-radius
positions; a 2e5 N/m attachment integrates with ordinary tolerances.

### Limits

No per-body aerodynamics, SRP or thermal loads, no contact, no attachment-to-attachment joints,
and no attachment on a spacecraft with the robot-arm effector. The `arm_*` robot-arm path is separate
and unchanged.

# Extensibility

Use this page when you are adding new models or control hooks to SpaceAGORA.

This page is for extension authors who want to stay on the stable root package
surface instead of binding new work to package implementation details.

Shortest successful starting point:

```text
templates/force_torque_model_template.jl
```

What to read next:

- [Public API](generated/public_api.md)
- [Concepts](user/concepts.md)
- [Maintainer Overview](maintainer/index.md)

SpaceAGORA now exposes a minimal stable extension contract from the root package.

The rule is simple:

1. subtype a stable abstract interface exported by `SpaceAGORA`
2. implement the matching root hook function
3. wire the model into the typed runtime configuration

Do not build new extensions by importing internal modules directly unless you
are intentionally working on SpaceAGORA internals.

## Add a force or torque model

Stable interface:

- `SpaceAGORA.AbstractForceTorqueModel`
- `SpaceAGORA.wrench`
- `SpaceAGORA.environment_requirements`
- `SpaceAGORA.solver_partition`
- `SpaceAGORA.gravity_backbone_structure`
- `SpaceAGORA.gravity_backbone_acceleration_ii`
- `SpaceAGORA.calcForceTorque`

Preferred additive methods:

```julia
SpaceAGORA.environment_requirements(model) -> SpaceAGORA.EffectorEnvironmentRequirements
SpaceAGORA.solver_partition(model) -> :explicit | :implicit
SpaceAGORA.gravity_backbone_structure(model) -> :unsupported | :position_only_static_gravity
SpaceAGORA.gravity_backbone_acceleration_ii(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> accel_ii
SpaceAGORA.wrench(model, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) -> (force_ii, torque_body)
```

Compatibility method:

```julia
SpaceAGORA.calcForceTorque(model, x, p, i) -> (force_n, torque_n_m)
```

Use `wrench` for new work when possible. The engine owns stage-consistent
sampling and passes typed state/environment data into the effector. Keep
`calcForceTorque` only for compatibility or when migrating existing models.
Override `solver_partition` only if the effector should move to the implicit
side of `split_imex`; the default is `:explicit`.
Override the gravity-backbone hooks only if the effector should participate in
`gravity_backbone_split`, either as a strict position-only static-gravity core
term or as an explicit translational velocity kick.

Solver partition contract:

- `split_imex` is the atmosphere-implicit IMEX path
- `multirate` remains the control-focused split path
- `gravity_backbone_split` is a fixed-step gravity-backbone split mode with a symplectic gravity core plus explicit SRP / N-body velocity kicks
- `gravity_backbone_split` remains translational-only; it is not a fully symplectic whole-system solve
- force returned by `wrench` is inertial-frame and torque is body-frame
- the current heat-load state is an accumulated heat-rate integral and stays explicit

Registration:

- add the model instance to `dynamic_effectors` in the simulation configuration

Template:

- `templates/force_torque_model_template.jl`

## Add a density model

Stable interface:

- `SpaceAGORA.AbstractDensityModel`
- `SpaceAGORA.getDensity`
- `SpaceAGORA.getDensityBatch!`

Required scalar method:

```julia
SpaceAGORA.getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)
```

Optional batch method:

```julia
SpaceAGORA.getDensityBatch!(rhos, Ts, winds, model, hs, lats, lons, el_time, wind, p)
```

Registration:

- set the model as `environment_model.density_model`

Template:

- `templates/density_model_template.jl`

## Add a control effector or hook

Stable interface:

- `SpaceAGORA.AbstractControlEffectorModel`
- `SpaceAGORA.calcControlEffect!`
- `SpaceAGORA.calcControlForceTorque`
- `SpaceAGORA.calcControlMassFlowRate`

Typical methods:

```julia
SpaceAGORA.calcControlEffect!(model, u, p, t, i)
SpaceAGORA.calcControlForceTorque(model, u, p, i, t)
SpaceAGORA.calcControlMassFlowRate(model, u, p, i, t)
```

Registration:

- add the model instance to `ControlModel(control_effectors=(...))`

Template:

- `templates/control_hook_template.jl`

## Planet and ephemerides extensions

Stable interfaces:

- `SpaceAGORA.AbstractPlanet`
- `SpaceAGORA.AbstractEphemeridesModel`

These interfaces are public because no-GRAM onboarding and future user-defined
body/back-end work need a stable contract, but the practical extension recipes
for custom planets and ephemerides are still thinner than the force/density/
control paths above. Treat them as advanced extension points.

## Design guidance

When adding a new model:

1. keep state in the model object, not in global variables
2. make units explicit in field names or docstrings
3. keep expensive caches out of hot dispatch paths unless they are typed
4. add a small typed constructor and one smoke test
5. document the model through the root `SpaceAGORA` hook surface, not by telling
   users to import internal modules

## Repository-owned examples

The templates above are intentionally small and live in the repository root:

- `templates/force_torque_model_template.jl`
- `templates/density_model_template.jl`
- `templates/control_hook_template.jl`

Use them as starting points, then add subsystem-specific tests close to the code
that owns the new model.

# Environment models and field kernels

The environment is what a simulation samples around the spacecraft: planet
properties, ephemerides, the atmosphere, terrain and the gravity field. This
folder supplies those values. The models that turn them into forces and torques
on a spacecraft are in `../dynamics/effectors/`, and so are the radiation and
magnetic-field models, which are written as effectors.

| Folder or file | Contents |
|---|---|
| `ephemerides/` | Planets and their shapes (`Planets`) and the ephemeris models, including the analytic `SimpleEphemeridesModel` (`EphemeridesModels`). |
| `physical_models.jl`, `atmosphere/` | `EnvironmentModels`: the density models, the frozen GRAM grid and near-surface models, the combined atmosphere, the named surrogate presets and the density epoch. |
| `terrain/` | Terrain models (`TerrainModels`). |
| `gravity/` | Analytic gravity-field kernels: central acceleration, Encke differencing and J2, J2 secular rates and the gravity gradient. |
| `gravity/gravity_effectors.jl`, `aerodynamics/aerodynamic_effectors.jl` | Small forwarding modules that re-export a few dynamic-effector names under their older paths. They hold no force-model code. |

The gravity kernels and the spacecraft side of the gravity models live in two
folders but one module. `gravity/gravity_models.jl` includes the kernels here
and the model types, run-time helpers and wrench adapters from
`../dynamics/effectors/gravity/` in a fixed order, all evaluated in
`DynamicEffectors.GravityEffectors`. Keep that order when you add a file, so
types and methods keep their identities and the numerical expressions do not
change.

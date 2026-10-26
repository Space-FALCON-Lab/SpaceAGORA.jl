# Interactive visualization ownership

The viewer reads saved trajectories and a scene sidecar after a run. It is
opt-in and changes nothing in the propagation equations, the atmosphere model
or the controllers.

| Owner | Responsibility |
| --- | --- |
| `src/interfaces/analysis/visualization/scene/` | Scene metadata, geometry conversion, body rotation samples, atmosphere overlays, saved-row export and the versioned viewer JSON |
| `src/models/vehicle/structure/mesh_geometry.jl` | Shared mesh reading, bounds, articulation and point sampling |
| `src/orchestration/simulation/callbacks/save_fields.jl` | Optional viewer fields saved alongside the existing results |
| `src/orchestration/simulation/engine/persistence.jl` | Writing the sidecar after a run that asked for it |
| `src/orchestration/simulation/campaigns/monte_carlo_visualization.jl` | Per-sample output directories and ensemble export over the existing campaign runner |
| `SpaceAGORAVisualizer.jl` (separate package) | HTML assembly, browser rendering, playback and standalone CSV/JSON import |
| `ext/SpaceAGORAVisualizerExt.jl` | The optional, versioned connection to the renderer; it owns no propagation |

`SimulationSettings.save_visualization_scene` (default `false`) turns the
viewer output on. With it off, the default save fields and the run behave as
before. With it on, `visualization_save_fields` adds the fields the viewer uses
when the run can supply them: thruster levels and directions, the Sun
direction, link and robot-arm poses, and density. If a run passes an explicit
`save_fields` list, the engine appends any of these the list lacks. Exporting a
bundle that is already saved does not run a simulation.

Several capabilities first built alongside the viewer now have owners of their
own and are outside this scene boundary: the GRAM epoch alignment hook
(`src/models/environment/atmosphere/density_epoch.jl`), the aerodynamic mesh
surrogate (`src/models/dynamics/effectors/aerodynamic_mesh_surrogate.jl`), the
Apollo-style descent guidance and control (`src/models/gnc/guidance/landing/`,
`src/models/gnc/control/landing/`), touchdown handling in the simulation
callbacks, plume-surface physics
(`src/models/dynamics/effectors/force_torque_models/plume_surface_interaction.jl`
with `src/orchestration/simulation/callbacks/plume_callbacks.jl`) and regional
terrain (`src/models/environment/terrain/`). The viewer only reads what those
models save.

The design is adapted from Evan Yu's original viewer contribution. For usage,
see the [visualization guide](../src/user/visualization.md); the
[SpaceAGORAVisualizer.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl)
repository documents the renderer and carries its third-party attributions.

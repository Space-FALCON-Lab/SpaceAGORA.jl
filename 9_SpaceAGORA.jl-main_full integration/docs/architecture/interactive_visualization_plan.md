# Interactive visualization ownership

The viewer reads saved trajectories and a scene sidecar. It is opt-in and
does not change the propagation equations, atmosphere model or controller.

| Owner | Responsibility |
| --- | --- |
| `src/interfaces/analysis/visualization/scene/` | Scene metadata, geometry conversion, body rotation samples, atmosphere overlays, saved-row export and versioned viewer JSON |
| `src/models/vehicle/structure/mesh_geometry.jl` | Shared mesh reading, bounds, articulation and point sampling |
| `src/orchestration/simulation/callbacks/save_fields.jl` | Optional viewer records sampled alongside the existing results |
| `src/orchestration/simulation/engine/persistence.jl` | Write the sidecar after a flagged run |
| `src/orchestration/simulation/campaigns/monte_carlo_visualization.jl` | Per-sample output directories and ensemble export over the existing campaign runner |
| `SpaceAGORAVisualizer.jl` | HTML assembly, browser rendering, playback and standalone CSV/JSON import |
| `ext/SpaceAGORAVisualizerExt.jl` | Optional, versioned renderer connection; no propagation ownership |

With visualization disabled, existing default save fields and run behavior
remain unchanged. Enabling it adds pose, density, Sun and reported thruster
fields when available. Explicit save-field lists receive missing visual
fields once. Exporting an already saved bundle does not run a simulation.

The features developed alongside the viewer in PR121 have since landed with
owners of their own and are not part of this scene boundary: the GRAM epoch
alignment hook (`src/models/environment/atmosphere/density_epoch.jl`), the
aerodynamic mesh surrogate (`src/models/dynamics/coupled/aerodynamic_mesh_surrogate.jl`),
the Apollo-style descent guidance and control
(`src/models/gnc/guidance/landing/`, `src/models/gnc/control/landing/`), touchdown handling
in the simulation callbacks, plume-surface physics
(`src/models/dynamics/coupled/force_torque_models/plume_surface_interaction.jl` with
`src/orchestration/simulation/callbacks/plume_callbacks.jl`) and regional terrain
(`src/models/environment/terrain/`). The viewer only reads what those models save.

Adapted from Evan Yu's PR121 at `c03d1dce`; see the user guide and
the independent SpaceAGORAVisualizer renderer guide for usage and third-party renderer attribution.
